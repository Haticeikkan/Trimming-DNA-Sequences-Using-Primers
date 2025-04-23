import os
from Bio import SeqIO
from Bio.Seq import Seq

# Primer setleri ve hedef bölgeler (ITS, Betatubulin, Calmodulin)
primers = {
    "ITS": {
        "forward": "TCCGTAGGTGAACCTGCGG",
        "reverse": "TCCTCCGCTTATTGATATGC"
    },
    "Betatubulin": {
        "forward": "GGTAACCAAATCGGTGCTGCTTTC",
        "reverse": "ACCCTCAGTGTAGTGACCCTTGGC"
    },
    "Calmodulin": {
        "forward": "CCGAGTACAAGGAGGCCTTC",
        "reverse": "CCGATAGAGGTCATAACGTGG"
    }
}

# Reverse primer dizilerini ters tamamlayıcıya çevir (5'→3' yönüne göre eşleşme yapılabilmesi için)
for region in primers:
    primers[region]["reverse"] = str(Seq(primers[region]["reverse"]).reverse_complement())

# Hamming mesafesi hesaplayan fonksiyon (eşleşmelerdeki uyuşmazlık sayısını bulur)
def hamming_distance(s1, s2):
    mismatches = 0
    for c1, c2 in zip(s1, s2):
        if c1 != c2 and c1 != 'N' and c2 != 'N':  # N bazları belirsiz olduğu için sayılmaz
            mismatches += 1
    return mismatches

# Verilen dizide, primer ile uyuşan tüm konumları (belirli toleransla) bulan fonksiyon
def find_with_mismatches_all(sequence, primer, max_mismatches):
    primer_len = len(primer)
    positions = []
    for i in range(len(sequence) - primer_len + 1):
        window = sequence[i:i+primer_len]
        mismatches = hamming_distance(window, primer)
        if mismatches <= max_mismatches:
            positions.append((i, mismatches))  # Eşleşme pozisyonu ve uyuşmazlık sayısı
    return positions

# Kullanıcıdan analiz yapılacak barcode klasörünü al
barcode_num = input("Hangi barcode klasörü analiz edilecek? (örn: barcode01): ").strip()
barcode_path = os.path.join("D:/ALSU/Trimming-DNA-Sequences-Using-Primers/RawData", barcode_num)

# Geçerli klasör var mı kontrol et
if not os.path.exists(barcode_path):
    print(f"❌ {barcode_path} klasörü bulunamadı.")
    exit()

# Klasördeki .fastq, .fasta veya .fa uzantılı dosyaları bul
valid_extensions = [".fastq", ".fasta", ".fa"]
sequence_files = [f for f in os.listdir(barcode_path) if any(f.endswith(ext) for ext in valid_extensions)]
if not sequence_files:
    print("❌ FASTQ veya FASTA dosyası bulunamadı.")
    exit()

# Her dosyayı sırayla işle
for filename in sequence_files:
    filepath = os.path.join(barcode_path, filename)
    file_ext = os.path.splitext(filename)[-1].lower()
    format_type = "fastq" if file_ext == ".fastq" else "fasta"

    try:
        records = list(SeqIO.parse(filepath, format_type))  # FASTA/FASTQ dosyasını oku
    except Exception as e:
        print(f"⚠️ {filename} okunurken hata oluştu: {e}")
        continue

    for record in records:
        sequence = str(record.seq)
        matched_regions = []  # Eşleşen bölgeleri tutmak için liste

        # Her bir hedef bölge için primer eşleşmelerini kontrol et
        for region, primer_pair in primers.items():
            fwd = primer_pair["forward"]
            rev = primer_pair["reverse"]

            max_fwd_mismatch = round(len(fwd) * 0.20)  # Maksimum %20 uyuşmazlık
            max_rev_mismatch = round(len(rev) * 0.20)

            # Senaryo 1: İleri primerden başlayıp ters primere kadar olan bölge
            fwd_matches = find_with_mismatches_all(sequence, fwd, max_fwd_mismatch)
            for fwd_index, fwd_mismatches in fwd_matches:
                search_region = sequence[fwd_index:]
                rev_matches = find_with_mismatches_all(search_region, rev, max_rev_mismatch)
                if rev_matches:
                    rev_index, rev_mismatches = rev_matches[0]
                    end_index = fwd_index + rev_index + len(rev)
                    subseq = sequence[fwd_index:end_index]
                    if 500 <= len(subseq) <= 700:  # Eşleşen dizinin uzunluğu kontrolü
                        matched_regions.append({
                            "region": region,
                            "subseq": subseq,
                            "fwd_mismatches": fwd_mismatches,
                            "rev_mismatches": rev_mismatches,
                            "fwd_mismatch_ratio": round(fwd_mismatches / len(fwd), 2),
                            "rev_mismatch_ratio": round(rev_mismatches / len(rev), 2)
                        })

            # Senaryo 2: Ters primerden başlayıp ileri primere kadar olan bölge
            rev_matches = find_with_mismatches_all(sequence, rev, max_rev_mismatch)
            for rev_index, rev_mismatches in rev_matches:
                search_region = sequence[rev_index:]
                fwd_matches = find_with_mismatches_all(search_region, fwd, max_fwd_mismatch)
                if fwd_matches:
                    fwd_index2, fwd_mismatches = fwd_matches[0]
                    end_index = rev_index + fwd_index2 + len(fwd)
                    subseq = sequence[rev_index:end_index]
                    if 500 <= len(subseq) <= 700:
                        matched_regions.append({
                            "region": region,
                            "subseq": subseq,
                            "fwd_mismatches": fwd_mismatches,
                            "rev_mismatches": rev_mismatches,
                            "fwd_mismatch_ratio": round(fwd_mismatches / len(fwd), 2),
                            "rev_mismatch_ratio": round(rev_mismatches / len(rev), 2)
                        })

        # Eşleşen bölgeleri FASTA dosyasına kaydet
        for idx, match in enumerate(matched_regions):
            output_dir = f"D:/ALSU/Trimming-DNA-Sequences-Using-Primers/{match['region']}_Output/{barcode_num}"
            os.makedirs(output_dir, exist_ok=True)

            out_filename = f"{record.id}_{match['region']}_{idx+1}.fasta"
            output_path = os.path.join(output_dir, out_filename)

            with open(output_path, "w") as out_f:
                out_f.write(f">{record.id}_{match['region']}_{idx+1}\n")
                for i in range(0, len(match['subseq']), 60):
                    out_f.write(match['subseq'][i:i+60] + "\n")  # 60 karakterlik satırlarla yaz

            # Kullanıcıya bilgi mesajı yazdır
            print(f"✅ {record.id} -> {match['region']} bölgesi kaydedildi → {out_filename}")
            print(f"🔸 Forward mismatches: {match['fwd_mismatches']} ({match['fwd_mismatch_ratio']*100:.0f}%)")
            print(f"🔸 Reverse mismatches: {match['rev_mismatches']} ({match['rev_mismatch_ratio']*100:.0f}%)")
