### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# ╔═╡ f9f9b4f2-62c9-11f0-3438-27cf7886b3aa
begin
	using Pkg
	# # Add from GitHub repository
	Pkg.add(url="https://github.com/symbiotic-engineering/MarineHydro.jl")
end

# ╔═╡ dea1bff4-7eb6-4204-acbf-c0140da4230d
begin
	#Intro to automatic differentiation 
	using ForwardDiff
	g(x) = sin(x) 
	
	x = 1.5
	println("g($x) = ", g(x))
	println("g'($x) = ", ForwardDiff.derivative(g, x))
	
end

# ╔═╡ cbc9e3ed-c8be-4713-9d04-d3e0ec9e007a
begin
	using MarineHydro
	using Zygote
	
end

# ╔═╡ a6fe44fb-c400-4207-b8c7-d243f735b9e7


# ╔═╡ 3919f706-fce7-4174-9e47-91e12b8b2c8c


# ╔═╡ 1ad76dfe-7fff-4856-a844-7a852e1020e3


# ╔═╡ b91b0dc0-da09-43e2-969c-be681e092c4f

function added_mass_program(radius,omega,dof)  
    mesh = differentiableMesh(radius) #fd
    A = calculate_radiation_forces(mesh,dof,omega)[1]
    return A
end

function damping_program(radius,omega,dof)  
    mesh = differentiableMesh(radius) #fd
    B = calculate_radiation_forces(mesh,dof,omega)[2]
    return B
end

function diffraction_program(ω,radius,dof)  
    mesh = differentiableMesh(radius) #fd
    force = DiffractionForce(mesh,ω,dof)
    return force
end

# ╔═╡ e004d41e-c4f1-4002-99df-1036f9a28dd4


# ╔═╡ Cell order:
# ╠═a6fe44fb-c400-4207-b8c7-d243f735b9e7
# ╠═dea1bff4-7eb6-4204-acbf-c0140da4230d
# ╠═3919f706-fce7-4174-9e47-91e12b8b2c8c
# ╠═f9f9b4f2-62c9-11f0-3438-27cf7886b3aa
# ╠═1ad76dfe-7fff-4856-a844-7a852e1020e3
# ╠═cbc9e3ed-c8be-4713-9d04-d3e0ec9e007a
# ╠═b91b0dc0-da09-43e2-969c-be681e092c4f
# ╠═e004d41e-c4f1-4002-99df-1036f9a28dd4
