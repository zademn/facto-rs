from Crypto.Util.number import getPrime
import numpy as np 
from random import getrandbits, randint

l = [16, 32, 64, 80]
for k in l:
    p = getPrime(k)
    q = getPrime(k)
    print(f'let p = Integer::from_str_radix("{p}", 10).unwrap();')
    print(f'let q = Integer::from_str_radix("{q}", 10).unwrap();')
    n = p * q
    print(f'let n = Integer::from_str_radix("{n}", 10).unwrap();')
    print()


