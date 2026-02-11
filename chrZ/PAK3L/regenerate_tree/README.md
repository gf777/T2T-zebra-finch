# Regenerate the tree as in Kong et al. 2010

Briefly, run these in succession (with the manuscript pdf as input) to get the NWK tree:

```
python extract_alignment_from_pdf.py input.pdf out.fa
python fix_alignment_and_write_formats.py out.fa pak_alignment_from_pdf.faa
python get_phy.py
python fix_names.py
```