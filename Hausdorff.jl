using LinearAlgebra, Distances, Random, StatsBase

"""Generates a random point in the box [a,b]."""
function RandomBox(a, b)
  n = length(a)
  return a + rand(n) .* (b - a)
end

"""Generates a random point in a Euclidean ball."""
function RandomBall(radius, center)
  n = length(center)
  x = randn(n)
  x = x / norm(x)
  r = rand()^(1 / n)
  return (radius * r) * x + center
end

"""Generates a random point in the intersection of a ball centered
at the origin and the nonnegative orthant."""
function RandomBallOrthant(radius, n)
  x = RandomBall(radius, zeros(n))
  return abs.(x)
end

"""Generates a random point in the probability simplex."""
function RandomSimplex(n)
  x = -log.(rand(n))  
  return x / sum(x)
end

"""Generates a random point in an L1 ball."""
function RandomL1Ball(radius, center)
  n = length(center)
  x = -log.(rand(n))  
  x = x / sum(x)
  for i = 1:n
    if rand() < 1 / 2
      x[i] = - x[i]
    end
  end
  r = rand()^(1 / n)
  return (radius * r) * x + center
end

"""Computes the Hausdorff distance between the point sets A and B."""
function hausdorff(A, B)
  D = pairwise(Euclidean(), A, B)
  dAB = maximum(minimum(D, dims = 2))
  dBA = maximum(minimum(D, dims = 1))
  return max(dAB, dBA)
end

"""Projects the point y onto a re-centered set."""
function RecenterProjection(Proj, y::Vector{T}, c::Vector{T}) where T <: Real
  return Proj(y - c) + c # set is translated by c
end

"""Projects the point y onto a scaled set."""
function ScaleProjection(Proj, y::Vector{T}, s::T) where T <: Real
  return s * Proj(y / s) # s > 0 is the scaling factor
end

"""Projects the point y onto the closed ball with radius r."""
function BallProjection(y::Vector{T}, r = one(T)) where T <: Real
#
  distance = norm(y)
  if distance > r
    return (r / distance) * y
  else
    return y
  end
end

"""Projects the point y onto the closed box with bounds a and b."""
function BoxProjection(y::Vector{T}, a = -ones(T, length(y)),
  b = ones(T, length(y))) where T <: Real
#
  return clamp.(y, a, b)
end

"""Projects the point y onto the simplex {x | x >= 0, sum(x) = r}."""
function SimplexProjection(y::Vector{T}, r = one(T)) where T <: Real
#
  n = length(y)
  z = sort(y, rev = true)
  (s, lambda) = (zero(T), zero(T))
  for i = 1:n
    s = s + z[i]
    lambda = (s - r) / i
    if i < n && lambda < z[i] && lambda >= z[i + 1]
      break
    end
  end
  return max.(y .- lambda, zero(T))
end

"""Projects the point y onto the ell_1 ball with radius r."""
function L1BallProjection(y::Vector{T}, r = one(T)) where T <: Real
#
  p = abs.(y)
  if norm(p, 1) <= r
    return y
  else
    x = SimplexProjection(p, r)
    return sign.(y) .* x
  end
end

"""Projects the point y onto the intersection of the ball of
radius r and the nonnegative orthant."""
function BallAndOrthantProjection(y::Vector{T}, r = one(T)) where T <: Real
#
  x = copy(y)
  x = max.(x, zero(T)) # project onto orthant
  return (r / max(norm(x), r)) .* x
end

"""Finds the support point for y on the inflated unit ball."""
function BallSupp(y::Vector{T}, r = one(T)) where T <: Real
#
  return (r / norm(y)) * y
end

"""Finds the support point for y on the box [a, b]."""
function BoxSupp(y::Vector{T}, a = -ones(T, length(y)),
  b = ones(T, length(y))) where T <: Real
#
  n = length(y)
  x = zeros(T, n)
  for i = 1:n
    if y[i] > zero(T)
      x[i] = b[i]
    elseif y[i] < zero(T)
      x[i] = a[i]
    else
      x[i] = (a[i] + b[i]) / 2
    end
  end
  return x
end

"""Finds the support point for y on the simplex {x | x >= 0, sum(x) = r}."""
function SimplexSupp(y::Vector{T}, r = one(T)) where T <: Real
#
  x = zeros(T, length(y))
  (v, m) = findmax(y)
  x[m] = r
  return x
end

"""Finds the support point for y on the L1 ball."""
function L1BallSupp(y::Vector{T}, r = one(T)) where T <: Real
#
  x = zeros(T, length(y))
  (v, m) = findmax(abs, y)
  x[m] = sign(y[m]) * r
  return x
end

"""Finds the support point for y on the intersection of the ball of
radius r and the nonnegative orthant."""
function BallAndOrthantSupp(y::Vector{T}, r = one(T)) where T <: Real
#
  x = max.(y, zero(T))
  if sum(x) <= zero(T)
    return zeros(T, length(y))
  else
    return (r / norm(x)) * x
  end
end

"""Projects the point y onto the Minkowski rounded set
R = c * S + (1 - c) * B. Here B is the unit ball, Proj 
is projection onto S, and Proj_R(y) = a + b."""
function MinkowskiNear(Proj, y, c, conv)
#
  n = length(y)
  (aold, bold) = (zeros(n), zeros(n))
  (anew, bnew) = (zeros(n), zeros(n))
  for iter = 1:100
    anew = c .* Proj((y - bold) ./ c) # project onto c * S
    bnew = (1 - c) .* BallProjection((y - anew) ./ (1 - c))
    if norm(aold - anew) + norm(bold - bnew) < conv
      break
    else
      @. aold = anew
      @. bold = bnew
    end
  end
  return anew + bnew
end

"""Finds the farthest point on A from B by Frank-Wolfe."""
function FrankWolfe(SuppA, PB, x0)
  (xold, xnew) = (copy(x0), similar(x0))
  for iter = 1:100
    xnew = SuppA(xold - PB(xold))
    if norm(xnew - xold) < 1.0e-10
      break
    else  
      xold .= xnew
    end
  end
  far = norm(xnew - PB(xnew))
  return (far, xnew)
end

"""Finds the farthest point on A from B by projected gradient ascent."""
function farthest(PA, PB, x0)
  (xold, xnew) = (copy(x0),copy(x0))
  for iter = 1:100
    xnew = PA(2xold - PB(xold))
    if norm(xnew - xold) < 1.0e-10 
      break
    else  
      xold .= xnew
    end
  end
  far = norm(xnew - PB(xnew))
  return (far, xnew)
end

"""Finds the farthest point on A from B by homotopy."""
function farthest_homotopy(PA, PB, SA, CenterA, CenterB, x0, n, method)
  x = BallProjection(x0)
  (far, homotopy_points, conv) = (0.0, 10, 1.0e-10) 
  for iter = 0:homotopy_points
    if iter == 0 # ball to ball
      d = norm(CenterA - CenterB)
      (far, x) = (d, (1 + 1 / d) * CenterA - CenterB / d)
    elseif iter == homotopy_points # d(A, B)
      if method == "proj grad" 
        (far, x) = farthest(PA, PB, x)
      elseif method == "Frank-Wolfe"
        (far, x) = FrankWolfe(SA, PB, x)
      end
    else # intermediate case
      t = iter / homotopy_points
      PMB(z) = MinkowskiNear(PB, z, t, conv)
      if method == "proj grad"
        PMA(z) = MinkowskiNear(PA, z, t, conv)
        (far, x) = farthest(PMA, PMB, x)
      elseif method == "Frank-Wolfe"
        SM(z) = SA(t * z) + BallSupp((1 - t) * z)
        (far, x) = FrankWolfe(SM, PMB, x)
      end
    end   
  end
  return (far, x) 
end

"""Orchestrates Hausdorff distance estimation."""
function master(ProjA, ProjB, SuppA, SuppB, CenterA, CenterB, method,
  homotopy, n, trials, io)
#
  (count, tries, optimum, obj) = (0, 100, 0.0, zeros(trials))
  x0 = zeros(n)
  PA(z) = RecenterProjection(ProjA, z, CenterA) 
  PB(z) = RecenterProjection(ProjB, z, CenterB)
  for trial = 1:trials
    success = false
    for i = 1:tries # find a point in A \ B
      x0 = PA(randn(n))
      if norm(PB(x0) - x0) > 1.0e-10
        success = true
        break
      end
    end
    if homotopy # solve for d(A, B)
      (objA, xA) = farthest_homotopy(PA, PB, SuppA, CenterA, CenterB, 
        x0, n, method)
    else
      if method == "proj grad" && success # solve for d(A, B)
        (objA, xA) = farthest(PA, PB, x0)
      elseif method == "Frank-Wolfe" && success
        (objA, xA) = FrankWolfe(SuppA, PB, x0)
      else
        objA = 0.0
      end
    end
    success = false
    for i = 1:tries # find point in B \ A
      x0 = PB(randn(n))
      if norm(PA(x0) - x0) > 1.0e-10
      success = true
        break
      end
    end
    if homotopy # solve for d(B, A)
      (objB, xB) = farthest_homotopy(PB, PA, SuppB, CenterB, CenterA,
         x0, n, method)
    else
      if method == "proj grad" && success
        (objB, xB,) = farthest(PB, PA, x0)
      elseif method == "Frank-Wolfe" && success
        (objB, xB) = FrankWolfe(SuppB, PA, x0)
      else
        objB = 0.0
      end
    end
    obj[trial] = max(objA, objB) # Hausdorff distance
    if obj[trial] > optimum + 10.0e-10 # update count of maximum distance
      count = 1
      optimum = obj[trial]
    elseif obj[trial] > optimum - 10.0e-8
      count = count + 1
    end
  end
  (avg, stdev) = (mean(obj), std(obj))
  if stdev < 1.0e-10 stdev = 0.0 end
  return (fraction = count / trials, optimum, avg, stdev)
end

outfile = "Hausdorff.out";
io = open(outfile, "w");
trials = 100;
points = 10000
println(io,"Set Pair"," & ","p"," & ","Method"," & ","Homotopy"," & ",
  "Maximum"," & ","Mean"," & ","Std"," &  ","Seconds"," \\ ")
for n in [2, 3, 10, 100, 1000]
  for i = 1:2
    if i == 1
      CenterA = ones(n)
      CenterB = zeros(n)
      ProjA = BoxProjection
      ProjB = BallAndOrthantProjection
      SuppA = BoxSupp
      SuppB = BallAndOrthantSupp
      title = "dH(box, ball and orthant)"
    elseif i == 2
      CenterA = zeros(n)
      CenterB = ones(n)
      ProjA = SimplexProjection
      ProjB = L1BallProjection
      SuppA = SimplexSupp
      SuppB = L1BallSupp
      title = "dH(simplex, L1 ball)"
    end
#
    (method, homotopy) = ("proj grad", false);
    Random.seed!(1234)
    time = @elapsed (fraction, optimum, avg, stdev) = master(ProjA, 
      ProjB, SuppA, SuppB, CenterA, CenterB, method, homotopy, n, 
      trials, io)
    println(io,title," & ",n," & ",method," & ",homotopy," & ",
        round(optimum, sigdigits=5)," & ",round(avg, sigdigits=5)," & ",
        round(stdev, sigdigits=5)," & ",round(time/trials, sigdigits=3)," \\ ")
#
    (method, homotopy) = ("proj grad", true);
    Random.seed!(1234)
    time = @elapsed (fraction, optimum, avg, stdev) = master(ProjA, 
      ProjB, SuppA, SuppB, CenterA, CenterB, method, homotopy, n, 
      trials, io)
    println(io,title," & ",n," & ",method," & ",homotopy," & ",
      round(optimum, sigdigits=5)," & ",round(avg, sigdigits=5)," & ",
      round(stdev, sigdigits=5)," & ",round(time/trials,sigdigits=3)," \\ ")
#
    (method, homotopy) = ("Frank-Wolfe", false);
    Random.seed!(1234)
    time = @elapsed (fraction, optimum, avg, stdev) = master(ProjA, 
      ProjB, SuppA, SuppB, CenterA, CenterB, method, homotopy, n, 
      trials, io)
    println(io,title," & ",n," & ",method," & ",homotopy," & ",
      round(optimum, sigdigits=5)," & ",round(avg, sigdigits=5)," & ",
      round(stdev, sigdigits=5)," & ",round(time/trials,sigdigits=3)," \\ ")
#
    (method, homotopy) = ("Frank-Wolfe", true);
    Random.seed!(1234)
    time = @elapsed (fraction, optimum, avg, stdev) = master(ProjA, 
      ProjB, SuppA, SuppB, CenterA, CenterB, method, homotopy, n, 
      trials, io)
    println(io,title," & ",n," & ",method," & ",homotopy," & ",
      round(optimum, sigdigits=5)," & ",round(avg, sigdigits=5)," & ",
      round(stdev, sigdigits=5)," & ",round(time/trials,sigdigits=3)," \\ ")
#
    (method, homotopy) = ("point cloud", false);
    Random.seed!(1234)
    points = 10000
    A = zeros(n, points)
    B = zeros(n, points)
    if i == 1
      (a, b) = (ones(n), 2 * ones(n))
      for j = 1:points
        A[:, j] = RandomBox(a, b)
        B[:, j] = RandomBallOrthant(1.0, n)
      end
    else
      for j = 1:points
        A[:, j] = RandomSimplex(n)
        B[:, j] = RandomL1Ball(1.0, ones(n))
      end
    end
    time = @elapsed optimum = hausdorff(A, B)
    println(io,title," & ",n," & ",method," & ",homotopy," & ",
      round(optimum, sigdigits=5)," & ",round(optimum, sigdigits=5)," & ",
      round(0.0, sigdigits=5)," & ",round(time, sigdigits=3)," \\ ")
  end
end
close(io)
