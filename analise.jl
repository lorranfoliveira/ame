module Analise

using LinearAlgebra

export No, Material, Elemento

# =======================================================================
"""
Implementa as propriedades de um nó.

id::Int64 -> Identificação do nó.
x::Float64 -> Coordenada x.
y::Float64 -> Coordenada y.
fx::Float64 -> Força na direção x.
fy::Float64 -> Força na direção y.
mxy::Float64 -> Momento fletor.
"""
struct No
    id::Int64
    x::Float64
    y::Float64
    fx::Float64
    fy::Float64
    mxy::Float64

    function No(id, x, y, fx=0, fy=0, mxy=0)
        if id < 1
            throw(ArgumentError("A identificação do nó deve ser mais que 1!"))
        end
        new(id, x, y, fx, fy, mxy)
    end
end

"""
Calcula a distância entre dois nós.

no₁::No -> Nó de referência 1.
no₂::No -> Nó de referência 2.
"""
function distancia(no₁::No, no₂::No)
    return norm([(no₂.x - no₁.x), (no₂.y - no₁.y)])
end

"""
Retorna os graus de liberdade do nó.
"""
function gl_no(no::No)::Array{Int64}
    return [3no.id - 2, 
            3no.id - 1, 
            3no.id]
end

"""
Retorna o vetor de forças nodais.
"""
function vetor_forcas(no::No)::Array{Float64}
    return [no.fx, no.fy, no.mxy]    
end

# =======================================================================
"""
Implementa as propriedades do material.

v: Coeficiente de Poisson.
E: Módulo de Elasticidade.
"""
struct Material
    id::Int64
    v::Float64
    E::Float64

    function Material(id, v, E)
        if id < 1
            throw(ArgumentError("O número de identificação do material não pode ser menor que 1!"))
        end

        if !(0 < v < 0.5)
            throw(ArgumentError("O Coeficiente de Poisson deve estar entre (0 < v < 0.5)!"))
        end

        if !(E > 0)
            throw(ArgumentError("O Módulo de Elasticidade deve ser maior que 0!"))
        end

        new(id, v, E)
    end
end

# =======================================================================
"""
Implementa as propriedades de uma elemento de pórtico.
"""
struct Elemento
    id::Int64
    no₁::No
    no₂::No
    area::Float64
    inercia::Float64
    material::Material

    function Elemento(id, no₁, no₂, area, inercia, material)
        if id < 1
            throw(ArgumentError("O número de identificação do elemento não pode ser menor que 1!"))
        end
        if !(area > 0)
            throw(ArgumentError("A área da seção deve ser maior que 0!"))
        end
        if !(inercia > 0)
            throw(ArgumentError("O momento de inércia da seção deve ser maior que 0!"))
        end
        new(id, no₁, no₂, area, inercia, material)
    end
end

"""
Comprimento do elemento.
"""
function comprimento(elemento::Elemento)
    return distancia(elemento.no₁, elemento.no₂)
end

"""
Retorna os graus de liberdade do elemento.
"""
function gl_elem(elemento::Elemento)::Array{Int64}
    gl_no₁ = gl_no(elemento.no₁)
    gl_no₂ = gl_no(elemento.no₂)
    return vcat(gl_no₁, gl_no₂)
end

"""
Ângulo de inclinação do elemento.
"""
function angulo(elemento::Elemento)
    dx = elemento.no₂.x - elemento.no₁.x
    dy = elemento.no₂.y - elemento.no₁.y

    return atan(dy / dx)
end

"""
Matriz de rotação do elemento.
"""
function mat_rot(θ)::Matrix{Float64}
    s = sin(θ)
    c = cos(θ)
    
    return [[c, s, 0, 0]
            [-s, c, 0, 0]
            [0, 0, c, s]
            [0, 0, -s, c]]
end

"""
Matriz de rigidez local do elemento.
"""
function ke(elemento::Elemento)
    L = comprimento(elemento)
    E = elemento.material.E
    I = elemento.inercia
    A = elemento.area

    k = [[E * A / L, 0, 0, -E * A / L, 0, 0]
    [0, 12E * I / L^3, 6E * I / L^2, 0, -12E * I / L^3, 6E * I / L^2]
    [0, 6E * I / L^2, 4E * I / L, 0, -6 * E * I / L^2, 2E * I / L]
    [-E * A / L, 0, 0, E * A / L, 0, 0]
    [0, -12E * I / L^3, -6E * I / L^2, 0, 12E * I / L^3, -6E * I / L^2]
    [0, 6E * I / L^2, 2E * I / L, 0, -6 * I / L^2,4E * I / L]]

    return k
end

"""
Matriz de rigidez do elemento rotacionada para o sistema global.
"""
function ke_rot(elemento::Elemento)
    k = ke(elemento)
    θ = angulo(elemento)
    r = mat_rot(θ)

    return r' * k * r
end

# =======================================================================
end

