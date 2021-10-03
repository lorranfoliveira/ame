include("analise.jl")

using Test
using .Analise

# =======================================================================
@testset "No" begin
    ### Fields ###
    no1 = No(1, 3, 2, 5, 6, 7)
    @test no1.id == 1
    @test no1.x == 3
    @test no1.y == 2
    @test no1.fx == 5
    @test no1.fy == 6
    @test no1.fxy == 7

    no2 = No(1, 3, 2)
    @test no2.fx == 0
    @test no2.fy == 0
    @test no2.fxy == 0

    ### Erros Fields ###
    @test_throws ArgumentError No(0, 0.2, 5e9)

    ### Métodos ###
    @test Analise.distancia(No(1, -30, -10), No(2, 40, 10)) ≈ 72.80109889

    @test Analise.gl_no(No(2, 0.2, 5e9)) == [4, 5, 6]

    @test Analise.vetor_forcas(No(1, 0.2, 5e9, 3, 4, 5.6)) == [3, 4, 5.6]

end

# =======================================================================
@testset "Material" begin
    # Fields
    m = Material(1, 5e9, 0.2)
    #id
    @test m.id == 1
    @test m.v == 0.2
    @test m.E == 5e9

    # id
    @test_throws ArgumentError Material(0, 1e9, 0.2)
    # E
    @test_throws ArgumentError Material(1, -1, 0.2)
    # v
    @test_throws ArgumentError Material(1, 1e9, -1)
    @test_throws ArgumentError Material(1, 1e9, 1)

    # Métodos
end

# =======================================================================
@testset "Elemento" begin
    no₁ = No(1, 3, 4, 2, 3, 6)
    no₂ = No(2, 5, 6, 5, 4, 7)
    mat = Material(1, 5e9, 0.2)

    ### Fields ###
    elem = Elemento(1, no₁, no₂, 10, 2, mat)

    @test elem.id == 1
    @test elem.no₁ ≡ no₁
    @test elem.no₂ ≡ no₂
    @test elem.area == 10
    @test elem.inercia == 2
    @test elem.material ≡ mat

    ### Erros Fields ###
    # id
    @test_throws ArgumentError Elemento(0, no₁, no₂, 10, 2, Material)
    # area
    @test_throws ArgumentError Elemento(1, no₁, no₂, 0, 2, Material)
    # inercia
    @test_throws ArgumentError Elemento(1, no₁, no₂, 10, 0, Material)

    ### Métodos ###
    @test Analise.comprimento(elem) ≈ 2.8284271247461903

    @test Analise.gl_elem(elem) == [1, 2, 3, 4, 5, 6]

    @test Analise.vetor_forcas(elem) == [2, 3, 6, 5, 4, 7]

    @test Analise.angulo(elem) ≈ π / 4

    # cos(π/4) == sin(π/4) == t
    t = sin(π / 4)
    kr = [t t 0 0 0 0
          -t t 0 0 0 0
          0 0 1 0 0 0
          0 0 0 t t 0
          0 0 0 -t t 0
          0 0 0 0 0 1]
    @test Analise.mat_rot(elem) ≈ kr


end

@testset "Estrutura" begin
    # Unidades: kN, cm
    no₁ = No(1, 0, 0, 0, 0, 0, true, true, true)
    no₂ = No(2, 200, 200, 1, -1, 100)
    no₃ = No(3, 400, 100, 0, 0, 0, false, true, false)
    nos = [no₁, no₂, no₃]

    m = Material(1, 2e5)

    area = 13.796
    inercia = 8.1052e2
    elm₁ = Elemento(1, no₁, no₂, area, inercia, m)
    elm₂ = Elemento(2, no₂, no₃, area, inercia, m)
    elementos = [elm₁, elm₂]

    estrutura = Estrutura(nos, elementos)

    ### Fields ###
    @test estrutura.nos ≡ nos
    @test estrutura.elementos ≡ elementos

    ### Métodos ###
    @test Analise.num_gls(estrutura) == 9
    
    @test Analise.vetor_forcas(estrutura) == [0, 0, 0, 1, -1, 100, 0, 0, 0]
    
    @test Analise.vetor_apoios(estrutura) == [1, 2, 3, 8]

    # t = Analise.ke_estrutura(estrutura)
    t = Analise.deslocamentos(estrutura)

end
