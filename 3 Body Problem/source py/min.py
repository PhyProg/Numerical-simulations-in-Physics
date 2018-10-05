d = []
fi = []
t = []
col = []
dd = 0
ff = 0
tt = 0
cc = 0
n = 0
N = 40000
d_min = 0.1

#scan

#N, d_min = map(float, input().split())

#print(N, d_min)
N = int(N)
s = input()
while (True):
    ff,dd,tt,cc = map(float, s.split())
    d.append(dd)
    fi.append(ff)
    t.append(tt)
    col.append(cc)
    n += 1
    try:
        s = input()
    except EOFError:
        break

#print(d, fi, t, col)

#print("Rezultati")

for i in range(0, n):
    if col[i] > 0: continue
    if i == 0:
        if d[i] <= d_min and d[i] < d[i+1] and abs(fi[i]-fi[i+1]) < 0.001:
            print(str(round(fi[i],10)) + " " + str(d[i]) + " " + str(t[i]) + " " + str(col[i]))
        continue

    if i == n-1:
        if d[i] <= d_min and d[i] < d[i-1] and abs(fi[i]-fi[i-1]) < 0.001:
            print(str(round(fi[i],10)) + " " + str(d[i]) + " " + str(t[i]) + " " + str(col[i]))
        continue

    if d[i] <= d_min and d[i] < d[i-1] and d[i] < d[i+1] and abs(fi[i]-fi[i+1]) < 0.001 and abs(fi[i]-fi[i-1]) < 0.001:
        print(str(round(fi[i],110)) + " " + str(d[i]) + " " + str(t[i]) + " " + str(col[i]))
