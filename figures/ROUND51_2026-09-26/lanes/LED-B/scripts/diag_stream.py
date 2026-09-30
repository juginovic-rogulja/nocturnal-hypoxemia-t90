import fitz, sys, re
d = fitz.open(sys.argv[1]); p = d[0]; cs = b"".join(d.xref_stream(x) for x in p.get_contents()).decode("latin-1")
i = cs.find("lineJoin"); print("first 'lineJoin' at", i, "context:", cs[max(0, i - 200):i + 60].replace("\n", " | "))
