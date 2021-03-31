#=
# lunar swing-by
=#

const dt = 1e-1
const num_of_output = 1000
const time_end = 5.0e6
const G = 6.672e-11

const day = 24.0 * 3600
const dtout = round(Int, time_end / (dt * num_of_output))

const mass_earth = 5.974e24
const mass_moon = 7.340e22
const orb_moon = 3.85e8
const velc_moon = sqrt(G * mass_earth / orb_moon)
const period_moon = 2 * pi * orb_moon * sqrt(orb_moon / (G * mass_earth))

const pos_launch = 6.371e6
const velc_launch = 1.114e4


mutable struct Object
    x::Float64
    y::Float64
    vx::Float64
    vy::Float64

    Object(x, y, vx, vy) = new(x, y, vx, vy)
    Object(x, y) = new(x, y, 0.0, 0.0)
    Object() = new(0.0, 0.0, 0.0, 0.0)
end

function move!(obj::Object, ax, ay)
    obj.vx += ax * dt
    obj.vy += ay * dt
    obj.x += obj.vx * dt
    obj.y += obj.vy * dt
end

function calc()
    data = Array{Float64, 2}(undef, 13, 0)
    time_launch = 0 * period_moon / 48
    moon_rad = (2 * pi / period_moon) * time_launch
    moon = Object(orb_moon * cos(moon_rad), orb_moon * sin(moon_rad), velc_moon * -sin(moon_rad), velc_moon * cos(moon_rad))
    sc = Object(pos_launch, 0.0, 0.0, velc_launch)
    time = 0.0
    etot_sc = 0.0

    n = 1
    while time < time_end
        r_me_sq = moon.x ^ 2 + moon.y ^ 2
        r_se_sq = sc.x ^ 2 + sc.y ^ 2
        r_sm_sq = (moon.x - sc.x) ^ 2 + (moon.y - sc.y) ^ 2
        r_me = sqrt(r_me_sq)
        r_se = sqrt(r_se_sq)
        r_sm = sqrt(r_sm_sq)

        a_me = -G * mass_earth / r_me_sq

        # moon
        move!(moon, a_me * moon.x / r_me, a_me * moon.y / r_me)

        # space craft
        a_c_x = -G * ((mass_earth * sc.x) / (r_se_sq * r_se) + (mass_moon * (sc.x - moon.x)) / (r_sm_sq * r_sm))
        a_c_y = -G * ((mass_earth * sc.y) / (r_se_sq * r_se) + (mass_moon * (sc.x - moon.y)) / (r_sm_sq * r_sm))
        move!(sc, a_c_x, a_c_y)

        # energy
        k = (sc.vx ^ 2 + sc.vy ^ 2) / 2
        u = -G * (mass_earth / r_se) - G * (mass_moon / r_sm)
        etot_sc = k + u

        time += dt
        n += 1

        if n % dtout == 1
            data = hcat(data, [n, time / day, moon.x, moon.y, moon.vx, moon.vy, sc.x, sc.y, sc.vx, sc.vy, r_se, r_sm, etot_sc])
        end
    end
    return data
end

function main()
    data = calc()
end

main()
