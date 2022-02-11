__all__ = ["fib"]
def fib(num):
    if num == 1:
        return 1
    elif num == 2:
        return 2
    else:
        return fib(num-1) + fib(num-2)

if __name__ == "__main__":
    seq = []
    num = 1
    while fib(num) <= 4000000:
        seq.append( fib(num) )
        num += 1
    even = [ n for n in seq if n % 2 == 0 ]
    print(sum(even))
