sum_squares = 0
sum_squared = 0
for n in range(1,100+1):
    sum_squares += n**2
    sum_squared += n
sum_squared *= sum_squared
print( sum_squared - sum_squares )
