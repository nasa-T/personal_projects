__all__ = ["factors", "primeFactors"]
def factors(num):
    return [ n for n in range(1,num+1) if num % n == 0 ]
def primeFactors(num):
    # f = factors(num)
    # return [ factors(n)[1] for n in f if len(factors(n)) == 2 ]
    n = num
    for pf in range( 2, int(n**(0.5)) + 1 ):
        while n % pf == 0:
            n /= pf
            if n == 1 or n == pf:
                return pf
if __name__ == "__main__":
    print( primeFactors(1851802432543) )
