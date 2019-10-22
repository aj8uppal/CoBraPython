def roundup(inp):
    outp = round(inp)
    return outp+int(outp < inp)
