using LinearAlgebra

"""
    gen(n, w, h, filename)

Writes an RBL file containing an alternating series of
circular arcs and other conic segments.
- `n` is the number of circular arcs
- `w` is the weight of the conic segments
- `h` is the height of the control point of the arcs
Example: gen(3, 0.4, 10, "test.rbl")
"""
function gen(n, w, h, filename)
    open(filename, "w") do f
        println(f, "$(2n) 2")
        for i in 1:n
            α = 2π * i / n
            p1 = [cos(α), sin(α), 0]
            α += π / n
            p2 = [cos(α), sin(α), 0]
            println(f, "$(p1[1]) $(p1[2]) $(p1[3]) 1")
            println(f, "0 0 0 $w")
            println(f, "$(p2[1]) $(p2[2]) $(p2[3]) 1")
            α += π / n
            p4 = [cos(α), sin(α), 0]
            p3 = (p2 + p4) / 2
            d = norm(p3 - p2)
            p3[3] = h
            println(f, "$(p2[1]) $(p2[2]) $(p2[3]) 1")
            println(f, "$(p3[1]) $(p3[2]) $(p3[3]) $(1 / √(h^2 / d^2 + 1))")
            println(f, "$(p4[1]) $(p4[2]) $(p4[3]) 1")
        end
    end
end
