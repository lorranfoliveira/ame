struct Quadrado
    a::Float64
    b::Float64

    function Quadrado(a::Float64, b::Float64)
        new(a, b)
    end
end

q = Quadrado(2, 3)