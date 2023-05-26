from falcon_gen.ntrugen import ntru_gen
from falcon_gen.ntt import div_zq

def falcon_gen(n):
  if n not in [ 2**i for i in range(1,11) ]:
    raise NotImplementedError("falconGen(n) supports only n = 2, 4, ..., 1024 but n = %d was given." % n)
  
  f, g, F, G = ntru_gen(n)
  h = div_zq(g, f)

  return f,g,h