def mult3_5(maximum):
    sum = 0
    for n in range(maximum):
        if (n % 3 == 0) or (n % 5 == 0):
            sum += n
    return sum