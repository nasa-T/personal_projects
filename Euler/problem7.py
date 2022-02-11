__all__ = ["isPrime"]

def isPrime(num):
    factors = []
    for n in range(1,num+1):
        if num % n == 0:
            factors.append(n)
        if len(factors) > 2:
            return False
    if len(factors) == 1:
        return False
    return True

if __name__ == "__main__":
    i = 0
    num = 1
    primes = []
    while True:
        if isPrime(num):
            i += 1
            primes.append(num)
            print(i)
        if i == 10002:
            print(primes[10001])
            break
        num += 1
