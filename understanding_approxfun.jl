f = Fun(Taylor, [1, 2, 3])
f(.1)

g(x) = 3 * x^2 + 2 * x + 1
g(.1)

f = Fun(Taylor() * Taylor(), [1, 2])
f(.1, .2)

g(x, y) = 2 * x^0 * y^1 + 1 * x^0 * y^0
g(.1, .2)

f = Fun(Taylor() * Taylor(), [1, 2, 3])
f(.1, .2)

g(x, y) =
    1 * x^0 * y^0 +
    2 * x^0 * y^1 +
    3 * x^1 * y^0

g(.1, .2)

coeforder(m) = [(j, d - j) for d in 0:m for j in 0:d]
coefindicesx(order, m = 5) = findall(t -> t[1] == order, coeforder(m))
coefindicesy(order, m = 5) = findall(t -> t[2] == order, coeforder(m))



binomial(9, 7)
index = findall(t -> t[1] == 2, coeforder(5))

geocoefs = out.out[6 : end-1]
coefdf = DataFrame(to = map(r -> sum(r), coeforder(7)),
                   coefs = geocoefs.array)
coefdf
avgs = combine(groupby(coefdf, :to), :coefs => mean)
scatter(coefdf.to, coefdf.coefs,
        xlab = "Total Order",
        ylab = "Geo coefs", label = "",
        title = "How do geo cheby coefs decline with total order")
plot!(avgs.to, avgs.coefs_mean, label = "")
