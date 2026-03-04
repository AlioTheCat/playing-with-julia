# Conservation laws :

# (1)  E₁ + E₂ = ½m₁(v₁)² + ½m₂(v₂)² = E ~ ½E
#
# E = m₁v²
# x = √m₁×v₁ / √E = v₁/v    ;   y = √m₂×v₂ / √E = (√m₂×v₂) / (√m₁×v) (but we don't really care)
#
#    ⟶  x² + y² = 1

# (2)  p₁ + p₂ = m₁v₁ + m₂v₂ = P
#
#    ⟶  √m₁(√E x) + √m₂(√E y) = P
#
#    ⟶  y = -(√m₁/√m₂)x + P/√E = -(√m₁/√m₂)x + P′

# (3) On collision, v₁ ⟵ -v₁

using Plots

global m₁ = BigFloat;
global m₂ = BigFloat;
global v = BigFloat;

global x = BigFloat;
global y = BigFloat;

function box_collision!(x::BigFloat, y::BigFloat)
    # Contraintes :
    # - y = αx + β
    # - x² + y² = 1

    α::BigFloat = -(√m₁ / √m₂)
    β = y - α*x

    # x² + y² = 1
    # Thereby :	(αx+β)² + x² = 1
    #			(α²+1).x² + 2αβ.x + (β²-1) = 0

    Δ = (2α*β)^2 - 4(α^2+1)*(β^2-1)
    @assert Δ ≥ 0

    # Modification des variables globales passées en argument
    global x = (-2α*β + √Δ) / ( 2(α^2 + 1) )
    global y = α*x + β

end

function wall_collision!(x::BigFloat, y::BigFloat)
    global y = -y
end

function can_catch_up(x,y)
    # the circle tangent line's slope -x/y should be greater
    # than the momentum conservation slope -√m₁/√m₂
    # so that the x coordinate can increase after a box collision
    # Also checks that Δ≥0, which is an equivalent criteria, but just for safety.
    α = -(√m₁ / √m₂)
    β = y - α*x

    Δ = (2α*β)^2 - 4(α^2+1)*(β^2-1)
    
    (x == -1) || ( 1-x > 0 && x/y ≤ √m₁/√m₂ && Δ ≥ 0)
end


function simulation(M=100, display=false, circle=true)
    global m₁ = M
    global m₂ = 1

    global x = BigFloat(-1) ;
    global y = BigFloat(0)

    # Now for the part we update :
    global N = 0;
    global X = Vector{Float64}([x]);
    global Y = Vector{Float64}([y]);

    while can_catch_up(x,y)

        #print(typeof(x), typeof(y))
        @assert typeof(x)==BigFloat
        @assert typeof(y)==BigFloat

        box_collision!(x,y)
        if display
            append!(X, Float64(x))
            append!(Y, Float64(y))
        end
        global N += 1

        if display
            print("Box collision ($N) : $x, $y \n")
        end

        if y < 0
            wall_collision!(x,y)
            if display
                append!(X, Float64(x))
                append!(Y, Float64(y))
            end
            global N += 1

            if display
                print("Wall collision ($N) : $x, $y \n")
            end
        end
    end
    print("The tiny box encountered $N collisions\n")
    
    if display && N ≤ 1000 # Let's keep this useful
        if circle
            fig = plot()
            plot!(fig, X, Y)
            t = range(0, 2π, 1000)
            plot!(fig, cos.(t), sin.(t)) 
            print("Plotted the states in the state space")
        else
            plot(X, Y)
        end
    end
end