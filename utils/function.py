
def s(x):
    if x <= -2 or x > 2:
        return 0
    elif x > -2 and x <= -1:
        return 1/4 * pow(2+x, 3)
    elif x > -1 and x <= 0:
        return 1/4 * (pow(2+x, 3) - 4 * pow(1+x, 3))
    elif x > 0 and x <= 1:
        return 1/4 * (pow(2-x, 3) - 4 * pow(1-x, 3))
    elif x > 1 and x <=2:
        return 1/4 * pow(2-x, 3)

def ds(x):
    if x <= -2 or x > 2:
        return 0
    elif x > -2 and x <= -1:
        return 1/4 * 3 * pow(2+x, 2)
    elif x > -1 and x <= 0:
        return 1/4 * (-12 * pow(1+x, 2) + 3 * pow(2+x, 2))
    elif x > 0 and x <= 1:
        return 1/4 * (12 * pow(1-x, 2) + 3 * pow(2-x, 2))
    elif x > 1 and x <=2:
        return 1/4 * -3 * pow(2-x, 2)

def fi(n,h,x,i):
    if i == 0:
        return s(x/h) - 4 * s((x+h)/h)
    elif i == 1:
        return s((x-h)/h) - s((x+h)/h)
    elif i >= 2 and i <= n-1:
        return s((x-i*h)/h)
    elif i == n:
        return s((x-n*h)/h) - s((x-(n+2)*h)/h)
    elif i == n+1:
        return s((x-(n+1)*h)/h) - 4 * s((x-(n+2)*h)/h)

def dfi(n, h, x, i):
   if i == 0:
       return ds(x / h) - 4 * ds((x + h) / h)
   elif i == 1:
       return ds((x - h) / h) - ds((x + h) / h)
   elif i >= 2 and i <= n - 1:
       return ds((x - i * h) / h)
   elif i == n:
        return ds((x - n * h) / h) - ds((x - (n + 2) * h) / h)
   elif i == n + 1:
       return ds((x - (n + 1) * h) / h) - 4 * ds((x - (n + 2) * h) / h)
