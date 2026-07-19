import pdfplumber, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
pdf = pdfplumber.open("Decoding Bispecific Antibody Developability Design Rules and Predictive Models from a 160-Member Library.pdf")
for i in range(8, min(20, len(pdf.pages))):
    text = pdf.pages[i].extract_text()
    if text:
        print(f"=== PAGE {i+1} ===")
        print(text)
        print()
