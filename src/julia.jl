using Unitful
using Unitful.DefaultSymbols
using Unitful: keV, h, c
using FFTW
FFTW.set_num_threads(Sys.CPU_THREADS)
using PyPlot

# a struct to hold the simulation info
struct PropagationInfo
    x
    y
    fx
    fy
end

# utilities
"The rectangle function."
rect(x) = one(x) .* (abs.(x) .< 0.5)

"frequencies for which the fft is calculated"
fftfreq(n) = ifftshift(-cld(n-1, 2) : fld(n-1, 2)) / n


# physics
"Transfer function for the Fresnel free propagation."
H(fx, fy, λ, z) = exp(im * 2π / λ * z) .* exp(-im * π * λ * z * (fx^2 + fy^2))
#H(fx, fy, λ, z) = exp(im * 2π / λ * z) .* exp(im * 2π/λ / z * (fx^2 + fy^2)) / (im * λ * z)

# gratings
singrating(x, y, m, L) = 0.5(1 + m * cos.(2π * x / L))
rectgrating(x, y, m, L; th = 0) = 0.5(1 + m * sign.(cos.(2π * x / L) - th))
rectgrating(pp::PropagationInfo, m, L; th = 0) = rectgrating.(
    pp.x, pp.y, m, L, th = th)
square(x, y, m, L) = m * rect(x / 2L) * rect(y / 2L)
circle(x, y, m, L) = m * rect(sqrt(x^2 + y^2) / 2L)
circle(pp::PropagationInfo, m, L) = circle.(pp.x, pp.y, m, L)
"Turns an amplitude grating into a phase one."
turn2phase(g) =  exp.(im * π / 2 * g)


function fresnelgrid(;δx = 0.1e-6, δy = δx,
        Nx = 2000, Ny = Nx, efficient = true)
    if efficient
        # FFTs are faster on arrays of lengths equal to powers of 2
        Nx = nextprod([2], Nx)
        Ny = nextprod([2], Ny)
    end

    # x are columns, y are rows
    x = δx * (collect(1:Nx) .- Nx/2)'
    y = δy * (collect(1:Ny) .- Ny/2)

    fx = (fftfreq(length(x)) ./ δx)'
    fy = fftfreq(length(y)) ./ δy

    return PropagationInfo(x, y, fx, fy)
end

tonormalspace(U) = ifft(U, (1,2))
tofourierspace(u) = fft(u, (1,2))
flatwave(pp) = ones(ComplexF64, length(pp.y), length(pp.x))

"Propagate a wavefront in the Fourier space"
propagate(U, fx, fy, λ, z) = U .* H.(fx, fy, λ, z)
propagate(U, pp::PropagationInfo, λ, z) = propagate(U, pp.fx, pp.fy, λ, z)

propagatenormal(u, pp::PropagationInfo, λ, z) = tonormalspace(
    propagate(tofourierspace(u), pp, λ, z))



pp = fresnelgrid(δx = 0.1μm, Nx = 4000, Ny = 1)
period = 10μm
g0 = rectgrating(pp, 1, period)
u1 = flatwave(pp) .* g0
λ = h * c / 46keV
u2(z) = propagatenormal(u1, pp, λ, z)


talbot_distance = 2period^2 / λ

z = reshape(collect(range(0mm, stop = 0.2talbot_distance, length = 2000)), 1, 1, :)
talbot_carpet = abs.(u2(z)).^2

plothunit = mm
plotvunit = μm

imshow(talbot_carpet[1, :, :],
    interpolation = "bicubic",
    origin="lower",
    extent=[
        ustrip(plothunit, z[1]), # left
        ustrip(plothunit, z[end]), # right
        ustrip(plotvunit, pp.x[1]), # bottom
        ustrip(plotvunit, pp.x[end]) # top
        ],
    aspect="auto")

xlabel("z position ($plothunit)")
ylabel("x position ($plotvunit)")

ylim(-20, 20)

savefig("$(dirname(dirname(@__FILE__)))/docs/carpets/julia.png")
