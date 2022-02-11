n = 20

while True:
    stop = True
    for d in range (1,20+1):
        if n % d != 0:
            stop = False
            break
    if stop:
        print(n)
        break
    n += 1
