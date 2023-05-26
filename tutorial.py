from lwe_with_hints import *

fileName = "tutorial/lwe_instance.json"
print("Loading LWE instance from", fileName)
A,b,q = loadLWEInstanceFromFile(fileName)

n,m = A.shape
print("The LWE parameters are n = %d, m = %d, q = %d." % (n,m,q))

print("Starting lattice reduction.")
lattice = LWELattice(A,b,q)
lattice.reduce()

print("Lattice reduction found the following vector:")
print(lattice.shortestVector)

print("The LWE secret is:")
print(lattice.s)

print("BKZ required the following blocksize:")
print(lattice.successBlocksize)

print("")
print("Restarting attack with hints.")

v_1 = vec( [1] + [0]*(n-1) )
v_2 = vec( [0] + [1] + [0]*(n-2) )
v_3 = vec( [0]*2 + [1] + [0]*(n-3) )
v_4 = vec( [0]*3 + [1] + [0]*(n-4) )
v_5 = vec( [-459, -441, 107, -207, 30, 358, -221, -483, 457, 96, 118, -241, 400, -478, 374, -46, -376, 415, 213, 476, -195, 25, -486, 444, 228, 313, -252, -182, -314, 105, -248, 163, 489, -388, 222, 110, -493, -491, 378, 213, 493, 48, 497, 138, 441, 140, 351, 135, -123, 414, -7, -344, -320, 54, 400, 230, -80, -85, -76, -475, 342, 276, 340, 1, 477, 158, -378, 146, 274, -355] )
v_6 = vec( [-315, 212, 432, 236, 423, -389, 67, -313, 365, 416, -180, -121, -472, 56, -468, -234, -305, 173, 444, -348, 261, -120, 249, -466, 491, -349, -140, 49, 99, 383, -321, 127, 447, -199, 443, -89, -384, -102, -106, 253, 495, 500, 8, -354, 115, 488, -88, 471, -291, 323, 349, 68, -253, -35, 34, -175, 172, -162, 123, -266, -6, 321, 481, -116, -60, 227, 163, -24, 56, 114] )
v_7 = vec( [1]*n )

lattice.integratePerfectHint( v_1, 0 )
lattice.integratePerfectHint( v_2, 0 )
lattice.integratePerfectHint( v_3, 2 )
lattice.integratePerfectHint( v_4, 1 )
lattice.integratePerfectHint( v_5, -1670 )
lattice.integratePerfectHint( v_6, 2381 )
lattice.integrateModularHint( v_7, 0, 2 )


print("Starting lattice reduction.")
lattice.reduce()

print("Lattice reduction found the following vector:")
print(lattice.shortestVector)

print("The LWE secret is:")
print(lattice.s)

print("BKZ required the following blocksize:")
print(lattice.successBlocksize)

