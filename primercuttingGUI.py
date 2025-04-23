# Gerekli kütüphaneler içe aktarılıyor
import os
import time
import threading
from tkinter import filedialog, messagebox
import customtkinter as ctk  # Gelişmiş bir tkinter teması
from Bio import SeqIO  # Biyolojik dizileri okumak için
from Bio.Seq import Seq  # Sekans işlemleri için

# Tema ayarları (arayüzün görünümü)
ctk.set_appearance_mode("light")
ctk.set_default_color_theme("blue")

# Kullanılacak primer setleri (her bir hedef bölge için ileri ve geri primerler)
primers = {
    "ITS": {
        "forward": "TCCGTAGGTGAACCTGCGG",
        "reverse": str(Seq("TCCTCCGCTTATTGATATGC").reverse_complement())  # Ters tamamlayıcısı alınır
    },
    "Betatubulin": {
        "forward": "GGTAACCAAATCGGTGCTGCTTTC",
        "reverse": str(Seq("ACCCTCAGTGTAGTGACCCTTGGC").reverse_complement())
    },
    "Calmodulin": {
        "forward": "CCGAGTACAAGGAGGCCTTC",
        "reverse": str(Seq("CCGATAGAGGTCATAACGTGG").reverse_complement())
    }
}

# Hamming mesafesi hesaplama fonksiyonu (kaç tane uyuşmayan karakter var)
def hamming_distance(s1, s2):
    mismatches = 0
    for c1, c2 in zip(s1, s2):
        if c1 != c2 and c1 != 'N' and c2 != 'N':  # 'N' belirsiz bazları temsil eder
            mismatches += 1
    return mismatches

# Primerin eşleştiği tüm pozisyonları belirleyen fonksiyon (belirli hata payıyla)
def find_with_mismatches_all(sequence, primer, max_mismatches):
    primer_len = len(primer)
    positions = []
    for i in range(len(sequence) - primer_len + 1):
        window = sequence[i:i+primer_len]
        mismatches = hamming_distance(window, primer)
        if mismatches <= max_mismatches:
            positions.append((i, mismatches))  # Eşleşme pozisyonu ve kaç hata olduğu
    return positions

# Ana analiz fonksiyonu
def analyze_fastq(barcode_path, mismatch_percentage, progress_var, log_box):
    # .fastq uzantılı dosyalar alınır
    fastq_files = [f for f in os.listdir(barcode_path) if f.endswith(".fastq")]
    total_files = len(fastq_files)
    start_time = time.time()

    for file_index, filename in enumerate(fastq_files):
        filepath = os.path.join(barcode_path, filename)
        try:
            records = list(SeqIO.parse(filepath, "fastq"))  # Dosya okunur
        except Exception as e:
            log_box.insert("end", f"⚠️ {filename} okunamadı: {e}\n")
            continue

        for record in records:
            sequence = str(record.seq)
            matched_regions = []

            for region, primer_pair in primers.items():
                fwd = primer_pair["forward"]
                rev = primer_pair["reverse"]

                # Maksimum izin verilen hata sayısı hesaplanır
                max_fwd_mismatch = round(len(fwd) * mismatch_percentage / 100)
                max_rev_mismatch = round(len(rev) * mismatch_percentage / 100)

                # 1. İleri → Geri tarama
                fwd_matches = find_with_mismatches_all(sequence, fwd, max_fwd_mismatch)
                for fwd_index, fwd_mismatches in fwd_matches:
                    search_region = sequence[fwd_index:]
                    rev_matches = find_with_mismatches_all(search_region, rev, max_rev_mismatch)
                    if rev_matches:
                        rev_index, rev_mismatches = rev_matches[0]
                        end_index = fwd_index + rev_index + len(rev)
                        subseq = sequence[fwd_index:end_index]
                        # Yalnızca 500-700 bp arası bölümler alınır
                        if 500 <= len(subseq) <= 700:
                            matched_regions.append((region, subseq, record.id, fwd_mismatches, rev_mismatches))

                # 2. Geri → İleri tarama
                rev_matches = find_with_mismatches_all(sequence, rev, max_rev_mismatch)
                for rev_index, rev_mismatches in rev_matches:
                    search_region = sequence[rev_index:]
                    fwd_matches = find_with_mismatches_all(search_region, fwd, max_fwd_mismatch)
                    if fwd_matches:
                        fwd_index2, fwd_mismatches = fwd_matches[0]
                        end_index = rev_index + fwd_index2 + len(fwd)
                        subseq = sequence[rev_index:end_index]
                        if 500 <= len(subseq) <= 700:
                            matched_regions.append((region, subseq, record.id, fwd_mismatches, rev_mismatches))

            # Elde edilen uygun bölümler kaydedilir
            for idx, (region, subseq, rec_id, fwd_mm, rev_mm) in enumerate(matched_regions):
                output_dir = os.path.join("D:/ALSU/Trimming-DNA-Sequences-Using-Primers", f"{region}_Output", os.path.basename(barcode_path))
                os.makedirs(output_dir, exist_ok=True)
                out_filename = f"{rec_id}_{region}_{idx+1}.fasta"
                output_path = os.path.join(output_dir, out_filename)

                with open(output_path, "w") as out_f:
                    out_f.write(f">{rec_id}_{region}_{idx+1}\n")
                    for i in range(0, len(subseq), 60):  # FASTA formatında 60 bp/satır
                        out_f.write(subseq[i:i+60] + "\n")

                log_box.insert("end", f"✅ {rec_id} → {region} kaydedildi: {out_filename}\n")

        # İlerleme ve kalan süre güncellenir
        progress = int(((file_index + 1) / total_files) * 100)
        progress_var.set(progress)
        elapsed = time.time() - start_time
        est_total = elapsed / (file_index + 1) * total_files
        remaining = est_total - elapsed
        log_box.insert("end", f"📊 {progress}% tamamlandı | Tahmini kalan süre: {int(remaining)} sn\n")

    messagebox.showinfo("Bitti", "✅ Tüm dosyalar başarıyla analiz edildi.")

# Analiz başlatma fonksiyonu (arayüzdeki butonla tetiklenir)
def start_analysis():
    path = filedialog.askdirectory(title="Barcode klasörünü seçin")
    if not path:
        return

    try:
        mismatch = float(mismatch_entry.get())
    except:
        messagebox.showerror("Hata", "Lütfen geçerli bir eşleşme yüzdesi girin")
        return

    log_box.delete("1.0", "end")  # Log kutusunu temizle
    thread = threading.Thread(target=analyze_fastq, args=(path, mismatch, progress_var, log_box))
    thread.start()  # Analizi arka planda başlat

# --- Kullanıcı Arayüzü (GUI) Tasarımı ---
app = ctk.CTk()
app.geometry("700x500")
app.title("DNA Bölge Analiz Aracı")

frame = ctk.CTkFrame(app, corner_radius=15)
frame.pack(padx=20, pady=20, fill="both", expand=True)

# Başlık etiketi
ctk.CTkLabel(frame, text="DNA Analiz Aracı", font=("Arial", 24)).pack(pady=(10, 20))

# Eşleşme yüzdesi giriş alanı
mismatch_label = ctk.CTkLabel(frame, text="Maksimum eşleşme yüzdesi (%):")
mismatch_label.pack()

mismatch_entry = ctk.CTkEntry(frame, placeholder_text="20")
mismatch_entry.insert(0, "20")
mismatch_entry.pack(pady=(0, 10))

# Başlat butonu
start_button = ctk.CTkButton(frame, text="Analizi Başlat", command=start_analysis)
start_button.pack(pady=(0, 10))

# İlerleme çubuğu
progress_var = ctk.IntVar()
progress_bar = ctk.CTkProgressBar(frame, variable=progress_var)
progress_bar.pack(fill="x", padx=20, pady=10)

# Log kutusu (sonuç ve hata mesajları burada gösterilir)
log_box = ctk.CTkTextbox(frame, height=200)
log_box.pack(padx=20, pady=10, fill="both", expand=True)

# Arayüzü başlat
app.mainloop()
