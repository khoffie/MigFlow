using ApproxFun, Plots

function chebytaylor(n)
    x = -1 : .01 : 1
    f(x) = 1 / (1 + 25x^2)
    taylor = Fun(f, Taylor(-1..1), n + 1)
    cheby = Fun(f, Chebyshev(-1..1), n + 1)
    p1 = plot(x, f.(x), label = "Runges function")
    plot!(x, real.(taylor.(x)), label = "taylor, degree $n",
          color = "red")
    plot!(x, real.(cheby.(x)), label = "cheby degree $n",
          color = "green")
    scatter!(points(taylor), 0.02ones(1), label = "Taylor points",
             color = "red")
    scatter!(points(cheby), -0.02ones(1), label = "Cheby points",
             color = "green")
    p2 = plot(x, f.(x) .- real.(taylor.(x)),
              title = "Error f(x) - g(x)",
              label = "Taylor degree $n error")
    plot!(x, f.(x) .- real.(cheby.(x)), label = "Cheby degree $n error")
    p3 = plot(x, (f.(x) .- real.(taylor.(x))) ./ f.(x),
              title = "Relative Error (f(x) - g(x)) / f(x)",
              label = "Taylor degree $n error")
    plot!(x, (f.(x) .- real.(cheby.(x))) ./ f.(x),
          label = "Cheby degree $n error")
    return p1, p2, p3
end


f(x) = 1 / (1 + 25x^2)
n = 15
taylor = Fun(f, Taylor(-1..1), n)
cheby = Fun(f, Chebyshev(-1..1), n)
p1 = plot(x, f.(x), label = "Runges function")
plot!(x, real.(taylor.(x)), label = "taylor, degree $n",
      color = "red")
plot!(x, real.(cheby.(x)), label = "cheby degree $n",
      color = "green")
scatter!(points(taylor), 0.02ones(1), label = "Taylor points",
         color = "red")
scatter!(points(cheby), -0.02ones(1), label = "Cheby points",
         color = "green")



p = chebytaylor(20)
p[3]
x = -1 : .01 : 1


using ApproxFun, Plots

p1 = plot(x, f.(x), label = "Runges function")
# degree = ncoefs - 1
plot!(x, real.(tay5.(x)), label = "degree 5", color = "red")
scatter!(points(tay5), zeros(6), label = "degree 5", color = "red")
plot!(x, real.(tay9.(x)), label = "degree 9", color = "green")
scatter!(points(tay9), zeros(10), label = "degree 9", color = "green")

points(cheby9)
points(tay5)

x = range(-2pi, 2pi, 100)
plot(x, cos.(x), label = "cos x", linewidth = 2)
plot!(x, cos.(2x), label = "cos(2x)")
plot!(x, cos.(0.5x), label = "cos(0.5x)")
