from pint import Quantity as Q_

a = Q_(1, 'm')
b = Q_(1, 's')

c = a / b
d = Q_(2, 'mph')
print((d + c).to('parsec/hour'))