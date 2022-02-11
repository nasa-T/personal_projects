def isPalindrome(s1):
    return s1 == s1[::-1]

palin_products = []
for n1 in range(100, 999 + 1):
    for n2 in range(100, 999 + 1):
        if isPalindrome(str(n1*n2)):
            palin_products.append(n1*n2)
print(max(palin_products))
