
using PyCall
using ForwardDiff
using LinearAlgebra
using MarineHydro
using ImplicitAD
using ChainRulesCore
include("meshGradients.jl")
cpt = pyimport("capytaine")

function bem(r1,dx1)
    rho = 1023
    omega= 1.03
    mesh = differentiableMeshPairs(r1,dx1)
    ptostiffness = -1e7
    C = hydrostatics(r1)
    M = hydrostatics(r1)
    Cmat = Diagonal([C,C])
    Mmat = Diagonal([M,M])
    total_nfaces = mesh.nfaces 
    # Define sphere classification functions - for linear solve and integration purposes
    # face_is_in_sphere_1(j) = j <= Int(total_nfaces/2)
    # face_is_in_sphere_2(j) = j > Int(total_nfaces/2)

    element_is_in_sphere_1(j) = get_center_x(j,mesh) <= r1 #assumes fixed resolution
    element_is_in_sphere_2(j) =  get_center_x(j,mesh) > r1 
    k = omega^2 / 9.81
    S, D = assemble_matrices((Rankine(), RankineReflected(), GFWu()), mesh, k)

    sphere_1_heave_normal = [element_is_in_sphere_1(j) ? mesh.normals[j,:]' * [0.0, 0.0, 1.0] : 0.0 for j in 1:total_nfaces]
    sphere_2_heave_normal = [element_is_in_sphere_2(j) ? mesh.normals[j,:]' * [0.0, 0.0, 1.0] : 0.0 for j in 1:total_nfaces]

    # sphere_1_center  = [element_is_in_sphere_1(j) ?  mesh.centers[j,:] : 0.0 for j in 1:total_nfaces ]
    # sphere_1_center = hcat(sphere_1_center'...)
    # sphere_2_center  = [element_is_in_sphere_2(j) ?  mesh.centers[j,:] : 0.0 for j in 1:total_nfaces]
    # sphere_2_center = hcat(sphere_2_center'...)

    # radiation of first sphere
    potential = MarineHydro.solve(D, S, -1im * omega * sphere_1_heave_normal)
    pressure = 1im * 1000 * omega * potential 
    force_on_sphere_1 = -sum(pressure .* sphere_1_heave_normal .* mesh.areas)
    A11 = real(force_on_sphere_1)/omega^2
    B11 = imag(force_on_sphere_1)/omega
    force_on_sphere_2 = -sum(pressure .* sphere_2_heave_normal .* mesh.areas)
    A12 = real(force_on_sphere_2)/omega^2
    B12 = imag(force_on_sphere_2)/omega
    Amat = [  #symmetric otherwise need bc with sphere_2_normal and itegrate on both sphere
        A11  A12;
        A12  A11
    ]
    Bmat = [
        B11  B12;
        B12  B11
    ]

    # @show Amat
    # @show Bmat
    inertia = Mmat+Amat  
    ptodamp = diag(Bmat)
   
  
  
    ##============Diffraction Force ====================== #
    dbc1 = -sum(airy_waves_velocity(mesh.centers,omega) .* mesh.normals, dims = 2) # diffraction boundary condition
    ϕd = solve(D,S,dbc1)
    dPressure = 1im*rho* ϕd * omega
    F1 = -sum(dPressure .* sphere_1_heave_normal .* mesh.areas)
    F2 = -sum(dPressure .* sphere_2_heave_normal .* mesh.areas)
    Airypressure1 = airy_waves_pressure(mesh.centers,  omega)
    Fk1 = -sum(Airypressure1 .* sphere_1_heave_normal .* mesh.areas)
    Airypressure2 = airy_waves_pressure(mesh.centers,  omega)
    Fk2 = -sum(Airypressure2 .* sphere_2_heave_normal .* mesh.areas)
    difraction_Force = [F1,F2]
    Froudekrylov = [Fk1,Fk2]
    return [A11, B11] #real(F1), imag(F1), real(Fk1), imag(Fk1), C, M #for one body for now
end

r1 = 1.0
dx1 = 4.0


@show bem(r1,dx1)
@show Zygote.jacobian(r1 -> bem(r1,dx1), r1)

