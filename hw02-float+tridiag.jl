### A Pluto.jl notebook ###
# v0.12.11

using Markdown
using InteractiveUtils

# ╔═╡ 60941eaa-1aea-11eb-1277-97b991548781
begin 
	using Pkg
	Pkg.activate(mktempdir())
	Pkg.add("PyPlot")
	Pkg.add("PlutoUI")
	Pkg.add("DataFrames")
	using PlutoUI
	using PyPlot
	using DataFrames
	using SparseArrays
	function pyplot(f;width=5,height=5)
		PyPlot.clf()
		PyPlot.grid()
		PyPlot.gcf().set_size_inches(width,height)
		PyPlot.gcf()
	end	
end

# ╔═╡ cec2ddf0-37f3-11eb-19d2-3791d3ebaf16
using LinearAlgebra

# ╔═╡ 4277c5e2-32f1-11eb-3bc9-73b0d499e912
md"""
# Homework II
Please return this assignment by  __Friday Dec. 11__. Fill in all the gaps 
marked with `...`  or `missing` and the empty plots. You may prototype things in another notebook or in the REPL, and cut & paste the results into this notebook.

Send the Pluto notebook and the corresponding pdf (created via the export button) in by e-mail to `juergen dot fuhrmann at wias-berlin dot de` or in a private message in the zulip chat. 


Please name the files according to your homework group, e.g. `HG-07-HW02.jl`,
`HG-07-HW02.pdf`


"""

# ╔═╡ e9f253c2-3577-11eb-2a9d-9d0688294a18
md"""
Homework group number: 12
"""

# ╔═╡ ff1cd722-3577-11eb-34ef-cb8f586801b0
md"""
Students: 
1) Aditya Kumar
2) Kay Töpfer
3) Tino Wagner
"""

# ╔═╡ d5425e00-32f1-11eb-0243-1fd9daf434fb
md"""
## 1) Summation (Basel sum)

We calculate $\sum_{n=1}^K \frac1{n^2}$ for different values of K and report the values for `Float16,Float32, Float64`:
"""

# ╔═╡ 3a07cec2-32f2-11eb-091e-a150b3d84726
exact_result=π^2/6

# ╔═╡ f8214cf6-32f1-11eb-2571-efc943d7fcf0
function simple_sum(T::Type,K)
	s=zero(T)
	for i=1:K
		ti=T(i)
		s=s+one(T)/(ti*ti)
	end
	s
end

# ╔═╡ 4ee29ce0-32fa-11eb-3e79-9769e01d5231
Ks=[10,100,1_000,10_000,100_000,1_000_000]

# ╔═╡ 9c4e40c4-32f2-11eb-143f-efee3b092e65
# Return a data frame which allows to report the data in a convenient way.
function basel_results(basel_sum::Function)
	result=DataFrame(K=Int[], Σ16=Float16[], error16=Float64[], Σ32=Float32[], error32=Float64[], Σ64=Float64[], error64=Float64[] )
	for K∈Ks
		Σ16=basel_sum(Float16,K)
		E16=abs(Σ16-exact_result)
		Σ32=basel_sum(Float32,K)
		E32=abs(Σ32-exact_result)
		Σ64=basel_sum(Float64,K)
		E64=abs(Σ64-exact_result)
		push!(result,(K,Σ16,E16,Σ32, E32,Σ64,E64 ))	
	end
	result
end

# ╔═╡ 8cdd12a8-32f2-11eb-07fb-afa1afa0212d
basel_results(simple_sum)

# ╔═╡ 749d952e-3580-11eb-13c4-c3fa59f3f6e8
md"""
What can be done to improve the accuracy of the calculation ?
"""

# ╔═╡ a706b49e-3b02-11eb-2bac-8757553291be
begin
	function more_accurate_sum(T::Type,K) 
		s = zero(BigFloat) 
		for k = 1:K
			k = BigInt(k)
		  	s += 1/k^2
		end
		convert(T,s)
		s
	end
end

# ╔═╡ adfbb620-3b02-11eb-3ccf-9b85c465eaa0
basel_results(more_accurate_sum)

# ╔═╡ 7f74d6f8-32f9-11eb-3c45-db4d4079e2fe
# Try to find a better way for summation here
function faster_sum(T::Type,K)
	s=collect(T,1:K).^-1
	s=dot(s,s)
end

# ╔═╡ b742f79a-32f9-11eb-17a1-ddd2d0b29a19
basel_results(faster_sum)

# ╔═╡ 16191448-32fa-11eb-1b43-b9bc0ee5ac91
# Write a report why your method is better
md"""
We found two solutions to increase accuracy. One is to simply increase K and the other would be to change the used data type T, so that rounding errors are decreased. We achieve this by using BigFloat. If we use BigFloat our run time increases.
The other solution (faster_sum) is  way faster (~10s versus ~500ms) and the error for Float32 and K=>10000 is slighly lesser. The use of vectorized calculation methods is much more preferable against loop involving methods.
"""

# ╔═╡ d8fe352e-357f-11eb-2335-2718b0fb72c7
md"""
__Bonus track :)__
"""

# ╔═╡ bd201d86-357f-11eb-186f-99e458460cf2
html"""
<iframe width="50%" height="200" src="https://www.youtube.com/embed/d-o3eB9sfls" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
"""

# ╔═╡ 491977ca-32f5-11eb-066e-23c1691688f7
md"""
## 2) 1D heat conduction problem

Given:

-  Domain $\Omega=(0,1)$
-  Right hand side $f: \Omega \to \mathbb{R}$, $f(x)=1$
-  Boundary values $v_L, v_R=0$
-  Transfer coefficient $\alpha=1$

Search function $u: \Omega \to \mathbb{R}$ such that

$-u'' =f \quad \text{in}\; \Omega$

$-u'(0) + \alpha (u(0)-v_L) = 0$

$u'(1) + \alpha (u(1)-v_R) = 1$
"""

# ╔═╡ 605cd878-3575-11eb-19a1-83daaa8ce42b
md"""
#### a) Calculate the exact solution of this problem:
"""

# ╔═╡ 75fb8be0-3575-11eb-15c3-6f2265dc2197
md"""
$-u''(x)=1$
$u(x) = -\frac{1}{2}x^2+c_1x+c_2$
$u'(x)= -x+c_1$
$-u'(0)+\alpha(u(0)-v_L)=0$
$c_1=c_2$
$u'(1)+\alpha(u(1)-v_R)=1$
$3c_1=\frac{5}{2}$
$c_1=\frac{5}{6}$
$c_2=\frac{5}{6}$
The exact solution is $u(x)= -\frac{1}{2}x^2+\frac{5}{6}x+\frac{5}{6}$
"""

# ╔═╡ 8e9a0f30-3575-11eb-237d-279b788492a6
md"""
#### b) What can you say about the limit of the solution for $\alpha\to\infty$?
"""

# ╔═╡ 09142fa2-3802-11eb-3967-ef51bfdf5a95
begin
	function c1(a)
		c2=(4+a)/(2a^2+4a)
		c1=a*c2
		c1
	end
	
	function c2(a)
		c2=(4+a)/(2a^2+4a)
		c2
	end
end

# ╔═╡ 3ff27110-3804-11eb-3830-6dd64275c7b6
let
X=collect(1:200)
PyPlot.clf() # Clear the figure
PyPlot.plot(X,map(c1,X),label="c1") # call the plot function
PyPlot.plot(X,map(c2,X),label="c2") # call the plot function
PyPlot.legend()
figure=PyPlot.gcf() # return figure to Pluto
end

# ╔═╡ a3fea48a-3575-11eb-22b8-0d995ac3d111
md"""
For $\alpha \rightarrow \infty$: $u(x)= -\frac{1}{2}x^2+\frac{1}{2}x$
"""

# ╔═╡ bfacdc24-3575-11eb-312e-37eeb3d3c192
md"""
#### c) Implement the finite difference discretization 

Use an equidistributed mesh with $N=2^{k}+1$:
"""

# ╔═╡ d065f9ba-3575-11eb-1895-77a505b8c0e4
k=10

# ╔═╡ ceac980e-3575-11eb-2b4c-17e953e1ea96
N=2^k+1

# ╔═╡ c860c236-3575-11eb-103e-997e6299c682
# Hint: you can use different functions for creating the vectors 
# of the tridiagonal matrix or a matrix made of these vectors
# If you create a full matrix here, you can extract diagonals via diag(A,i)
function heatmatrix(N,α)
	h=1/(N-1)
	A=Tridiagonal(fill(-1/h, N-1),fill(2/h, N),fill(-1/h, N-1))
	A[1,1]=1/h + α
	A[N,N]=1/h + α
	A
end

# ╔═╡ 152cc650-3576-11eb-2e8e-f1373219482c
function rhs(N,α)
	h=1/(N-1)
	f=fill(h, N)
	f[1]=h/2
	f[N]=h/2+1
	f
end

# ╔═╡ 03a9c0cc-3576-11eb-0a2f-23965446f5d0
md"""
#### d) Compare different solution strategies 
  
- Your implementation of TDMA (Progonka)
    - Julia dense matrix LU factorization
    - Julia tridiagonal matrix LU factorization
    - Julia sparse matrix LU factorization
    - Multiplication by the inverse

"""

# ╔═╡ 36be29f0-3830-11eb-3af2-b7b4dc2540c3
function solve_tdma(A,b)
	N=length(b)
	#forward sweep
	alpha=zeros(N-1)
	beta=zeros(N)
	alpha[1]=A[1,2]/A[1,1]
	beta[1]=b[1]/A[1,1]
	for i=2:N-1
		alpha[i]=A[i,i+1]/(A[i,i]-alpha[i-1]*A[i-1,i])
		beta[i]=(b[i]-A[i,i-1]*beta[i-1])/(A[i,i]-A[i,i-1]*alpha[i-1])
	end
	#backward sweep
	beta[N]=(b[N]-A[N,N-1]*beta[N-1])/(A[N,N]-A[N,N-1]*alpha[N-1])
	u=zeros(N)
	u[N]=beta[N] 
	for i=1:N-1
		i=N-i
		u[i]=beta[i]-alpha[i]*u[i+1]
	end
	u
end

# ╔═╡ 8983ffc0-383a-11eb-3d8a-b9bbe70cc0a0
function solve_julia_dense(A,b)
	N=length(b)
	A=zeros(N,N)+A
	lu(A)\b
end

# ╔═╡ 5c659588-3576-11eb-2c5a-4749ae27b97c
function solve_julia_triag(A,b)
	lu(A)\b
end

# ╔═╡ 6d297880-3576-11eb-0186-c5e21eb76e13
function solve_julia_sparse(A,b)
	A=sparse(A)
	lu(A)\b
end

# ╔═╡ 760538e0-3576-11eb-2e71-fd4b9cfdc299
function solve_inv_multiply(A,b)
	u=inv(A)*b
	u
end

# ╔═╡ bef8ba2c-3588-11eb-1a02-7b21323e43f8
md"""
Solve for different values of $k$, report the timings, e.g. using `@elapsed`
"""

# ╔═╡ 9ced48d0-3869-11eb-0a42-e34b1534bd22
begin
	function timings(ks)
result=DataFrame(k=Int[],timing_tdma=Float16[],timing_dense=Float16[],timing_triag=Float16[],timing_sparse=Float16[],timing_inv=Float16[] )
		for k∈ks
			N=2^k+1
			A=heatmatrix(N,1)
			b=rhs(N,1)
			timing_tdma=@elapsed solve_tdma(A,b)
			timing_dense=@elapsed solve_julia_dense(A,b)
			timing_triag=@elapsed solve_julia_triag(A,b)
			timing_sparse=@elapsed solve_julia_sparse(A,b)
			timing_inv=@elapsed solve_inv_multiply(A,b)
			push!(result,(k,timing_tdma,timing_dense,timing_triag,timing_sparse,timing_inv))
		end
		result
	end
end

# ╔═╡ d7dd8ef2-3588-11eb-07e5-c16aa053a910
md"""
Timings: see table down below. All values (except k) are in seconds.
"""

# ╔═╡ a4ed4620-3873-11eb-2c04-19dae0d6b212
timings(collect(1:12))

# ╔═╡ c3957852-3576-11eb-2c08-e11308b90ca7
md"""
#### e) Plot graphics of the solution for different resolution steps $k$:
"""

# ╔═╡ 1da13660-386c-11eb-0216-fbbbac3e339a
begin
	A=heatmatrix(N,1)
	b=rhs(N,1)
end

# ╔═╡ 95f600cc-3589-11eb-02c3-d5339e136f72
md"""
Discussion: Comparing the diffrent solution strategies, they all give an equivalent result (up to the scale of errors we were looking at). The approximation is getting smoother for higher k.
"""

# ╔═╡ 03c46ad4-3577-11eb-0c1c-87c0f1f83914
md"""
#### f) Plot the error etween exact and approximate solution vs. $h$

(Hint: use a log-log representation):
"""

# ╔═╡ 8d2cf9a0-3831-11eb-1398-913292dcc6f4
function exact_solution(x) #alpha = 1
	u=-(1/2)*x^2+(5/6)*x+(5/6)
	u
end

# ╔═╡ 55f2e0d8-3588-11eb-0d89-c32b17fa1385
begin
	function plot_results(k)
		N=2^k+1
		h=1/(N-1)
		x=collect(0:h:1)
		A=heatmatrix(N,1)
		b=rhs(N,1)
		PyPlot.clf() # Clear the figure
		PyPlot.plot(x,solve_tdma(A,b),label="tdma") # call the plot function
		PyPlot.plot(x,solve_julia_dense(A,b),label="dense") # call the plot function
		PyPlot.plot(x,solve_julia_triag(A,b),label="triag") # call the plot function
		PyPlot.plot(x,solve_julia_sparse(A,b),label="sparse") # call the plot function
		PyPlot.plot(x,solve_inv_multiply(A,b),label="multiply") # call the plot function
		PyPlot.plot(x,map(exact_solution, x),label="exact")
		PyPlot.legend()
		PyPlot.title("k = $k")
		figure=PyPlot.gcf() # return figure to Pluto
	end
end

# ╔═╡ 10f3c680-3876-11eb-3fa5-03b37a3cdffc
begin
	function get_errors(ks)
errors=DataFrame(h=Float16[],error_tdma=Float32[],error_dense=Float32[],error_triag=Float32[],error_sparse=Float32[],error_inv=Float32[] )
		for k=ks
			N=2^k+1
			h=1/(N-1)
			x=collect(0:h:1)
			A=heatmatrix(N,1)
			b=rhs(N,1)
			error_tdma=norm(map(exact_solution,x)-solve_tdma(A,b), Inf)
			error_dense=norm(map(exact_solution,x)-solve_julia_dense(A,b), Inf)
			error_triag=norm(map(exact_solution,x)-solve_julia_triag(A,b), Inf)
			error_sparse=norm(map(exact_solution,x)-solve_julia_sparse(A,b), Inf)
			error_inv=norm(map(exact_solution,x)-solve_inv_multiply(A,b), Inf)
			push!(errors,(h,error_tdma,error_dense,error_triag,error_sparse,error_inv))
		end
		errors
	end
	function plot_errors(ks)
		errors=Matrix(get_errors(ks))
		h=errors[:,1]
		PyPlot.clf() # Clear the figure
		PyPlot.loglog(h,errors[:,2],label="tdma") # call the plot function
		PyPlot.loglog(h,errors[:,3],label="dense") # call the plot function
		PyPlot.loglog(h,errors[:,4],label="triag") # call the plot function
		PyPlot.loglog(h,errors[:,5],label="sparse") # call the plot function
		PyPlot.loglog(h,errors[:,6],label="multiply") # call the plot function
		PyPlot.legend()
		PyPlot.ylabel("error")
		PyPlot.xlabel("h")
		figure=PyPlot.gcf() # return figure to Pluto
	end
end

# ╔═╡ dfa03fd0-3881-11eb-3390-d3c11fce6442
plot_errors(collect(1:12))

# ╔═╡ 9e93929e-3589-11eb-2238-233d3b611708
md"""
Discussion: For smaller grid size (h), the erros are getting higher for $\alpha = 1$. 
The finite difference method should increase accuracy by adding grid points. But, we noticed that with the increase of grid points, there is also an increase of rounding errors for the very small decimal values.
"""

# ╔═╡ c57b37ee-3576-11eb-2483-856e4faa572e
md"""
#### g) What happens for increased values of the transfer coefficient ?
  
Provise results e.g. for
$\alpha=1,10,100,1.0\cdot 10^5, 1.0\cdot 10^{10}, 1.0\cdot 10^{20}$
"""

# ╔═╡ a97f3844-3588-11eb-14c7-bf264804c4a1
begin
	function plot_results(k,alpha)
		N=2^k+1
		x=collect(0:1/(N-1):1)
		A=heatmatrix(N,alpha)
		b=rhs(N,alpha)
		PyPlot.clf() # Clear the figure
		PyPlot.plot(x,solve_tdma(A,b),label="tdma") # call the plot function
		PyPlot.plot(x,solve_julia_dense(A,b),label="dense") # call the plot function
		PyPlot.plot(x,solve_julia_triag(A,b),label="triag") # call the plot function
		PyPlot.plot(x,solve_julia_sparse(A,b),label="sparse") # call the plot function
		PyPlot.plot(x,solve_inv_multiply(A,b),label="multiply") # call the plot function
		PyPlot.legend()
		PyPlot.title("alpha = $alpha , k = $k")
		figure=PyPlot.gcf() # return figure to Pluto
	end
end

# ╔═╡ 27be7290-386f-11eb-0d89-f9d965d6f91f
plot_results(1)

# ╔═╡ 81323a50-386f-11eb-178e-b3c65b34a7a7
plot_results(2)

# ╔═╡ 2b4050e0-3870-11eb-2f2b-d50607778688
plot_results(5)

# ╔═╡ c4fd66d0-3881-11eb-2e7c-e1859d2465fb
plot_results(10)

# ╔═╡ 108544e0-3884-11eb-2495-4f3e0032446f
plot_results(10,1)

# ╔═╡ 1c397c20-3884-11eb-23be-a1dc06298fb5
plot_results(10,100)

# ╔═╡ 28b39b70-3884-11eb-0ee8-eb67d32ec9f4
plot_results(10,10^5)

# ╔═╡ 386fc660-3884-11eb-2a32-1ba13f7e789a
plot_results(10,10^20)

# ╔═╡ b6d3b92c-3589-11eb-15d6-57a482157251
md"""
Discussion: The solution approaches $-\frac{1}{2}x^2+\frac{1}{2}x$ for large $\alpha$, like we already concluded in 2b).
"""

# ╔═╡ Cell order:
# ╠═60941eaa-1aea-11eb-1277-97b991548781
# ╠═4277c5e2-32f1-11eb-3bc9-73b0d499e912
# ╠═e9f253c2-3577-11eb-2a9d-9d0688294a18
# ╠═ff1cd722-3577-11eb-34ef-cb8f586801b0
# ╟─d5425e00-32f1-11eb-0243-1fd9daf434fb
# ╠═3a07cec2-32f2-11eb-091e-a150b3d84726
# ╠═f8214cf6-32f1-11eb-2571-efc943d7fcf0
# ╠═4ee29ce0-32fa-11eb-3e79-9769e01d5231
# ╠═9c4e40c4-32f2-11eb-143f-efee3b092e65
# ╠═8cdd12a8-32f2-11eb-07fb-afa1afa0212d
# ╟─749d952e-3580-11eb-13c4-c3fa59f3f6e8
# ╠═cec2ddf0-37f3-11eb-19d2-3791d3ebaf16
# ╠═a706b49e-3b02-11eb-2bac-8757553291be
# ╠═adfbb620-3b02-11eb-3ccf-9b85c465eaa0
# ╠═7f74d6f8-32f9-11eb-3c45-db4d4079e2fe
# ╠═b742f79a-32f9-11eb-17a1-ddd2d0b29a19
# ╠═16191448-32fa-11eb-1b43-b9bc0ee5ac91
# ╟─d8fe352e-357f-11eb-2335-2718b0fb72c7
# ╟─bd201d86-357f-11eb-186f-99e458460cf2
# ╟─491977ca-32f5-11eb-066e-23c1691688f7
# ╟─605cd878-3575-11eb-19a1-83daaa8ce42b
# ╠═75fb8be0-3575-11eb-15c3-6f2265dc2197
# ╟─8e9a0f30-3575-11eb-237d-279b788492a6
# ╠═09142fa2-3802-11eb-3967-ef51bfdf5a95
# ╠═3ff27110-3804-11eb-3830-6dd64275c7b6
# ╠═a3fea48a-3575-11eb-22b8-0d995ac3d111
# ╟─bfacdc24-3575-11eb-312e-37eeb3d3c192
# ╠═d065f9ba-3575-11eb-1895-77a505b8c0e4
# ╠═ceac980e-3575-11eb-2b4c-17e953e1ea96
# ╠═c860c236-3575-11eb-103e-997e6299c682
# ╠═152cc650-3576-11eb-2e8e-f1373219482c
# ╟─03a9c0cc-3576-11eb-0a2f-23965446f5d0
# ╠═36be29f0-3830-11eb-3af2-b7b4dc2540c3
# ╠═8983ffc0-383a-11eb-3d8a-b9bbe70cc0a0
# ╠═5c659588-3576-11eb-2c5a-4749ae27b97c
# ╠═6d297880-3576-11eb-0186-c5e21eb76e13
# ╠═760538e0-3576-11eb-2e71-fd4b9cfdc299
# ╠═bef8ba2c-3588-11eb-1a02-7b21323e43f8
# ╠═9ced48d0-3869-11eb-0a42-e34b1534bd22
# ╠═d7dd8ef2-3588-11eb-07e5-c16aa053a910
# ╠═a4ed4620-3873-11eb-2c04-19dae0d6b212
# ╟─c3957852-3576-11eb-2c08-e11308b90ca7
# ╟─1da13660-386c-11eb-0216-fbbbac3e339a
# ╠═55f2e0d8-3588-11eb-0d89-c32b17fa1385
# ╠═27be7290-386f-11eb-0d89-f9d965d6f91f
# ╠═81323a50-386f-11eb-178e-b3c65b34a7a7
# ╠═2b4050e0-3870-11eb-2f2b-d50607778688
# ╠═c4fd66d0-3881-11eb-2e7c-e1859d2465fb
# ╠═95f600cc-3589-11eb-02c3-d5339e136f72
# ╟─03c46ad4-3577-11eb-0c1c-87c0f1f83914
# ╠═8d2cf9a0-3831-11eb-1398-913292dcc6f4
# ╠═10f3c680-3876-11eb-3fa5-03b37a3cdffc
# ╠═dfa03fd0-3881-11eb-3390-d3c11fce6442
# ╠═9e93929e-3589-11eb-2238-233d3b611708
# ╟─c57b37ee-3576-11eb-2483-856e4faa572e
# ╠═a97f3844-3588-11eb-14c7-bf264804c4a1
# ╠═108544e0-3884-11eb-2495-4f3e0032446f
# ╠═1c397c20-3884-11eb-23be-a1dc06298fb5
# ╠═28b39b70-3884-11eb-0ee8-eb67d32ec9f4
# ╠═386fc660-3884-11eb-2a32-1ba13f7e789a
# ╠═b6d3b92c-3589-11eb-15d6-57a482157251
