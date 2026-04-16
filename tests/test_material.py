import pytest
from src.materials.material import Material

def test_material_initialization():
    mat = Material(idt=1, young=210e9, poisson=0.3, density=7850.0)
    assert mat.idt == 1
    assert mat.young == 210e9
    assert mat.poisson == 0.3
    assert mat.density == 7850.0

def test_material_shear_modulus():
    # G = E / (2 * (1 + v))
    # Para E=210e9, v=0.3: G = 210e9 / 2.6 = 80.769e9
    mat = Material(idt=1, young=210e9, poisson=0.3, density=7850.0)
    assert mat.shear_modulus == pytest.approx(80.76923076923077e9)

def test_material_bulk_modulus():
    # K = E / (3 * (1 - 2v))
    # Para E=210e9, v=0.3: K = 210e9 / (3 * 0.4) = 210e9 / 1.2 = 175e9
    mat = Material(idt=1, young=210e9, poisson=0.3, density=7850.0)
    assert mat.bulk_modulus == pytest.approx(175e9)

def test_material_to_dict():
    mat = Material(idt=1, young=210e9, poisson=0.3, density=7850.0)
    expected = {
        "idt": 1,
        "young": 210e9,
        "poisson": 0.3,
        "density": 7850.0,
        "shear_modulus": mat.shear_modulus,
        "bulk_modulus": mat.bulk_modulus
    }
    assert mat.to_dict() == expected

def test_material_validation():
    with pytest.raises(TypeError, match="Young's modulus must be a number"):
        Material(idt=1, young="x", poisson=0.3, density=7850.0)
    
    with pytest.raises(ValueError, match="Young's modulus must be a positive number"):
        Material(idt=1, young=-100, poisson=0.3, density=7850.0)
        
    with pytest.raises(ValueError, match="Poisson's ratio must be between 0 and 0.5"):
        Material(idt=1, young=210e9, poisson=0.6, density=7850.0)

    with pytest.raises(ValueError, match="Density must be a positive number"):
        Material(idt=1, young=210e9, poisson=0.3, density=0)
