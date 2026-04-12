import pytest
import numpy as np

from src.fea.section import RectangularSection


def test_rectangular_section_initialization():
    section = RectangularSection(idt=1, b=2.0, h=3.0)
    assert section.idt == 1
    assert section.b == 2.0
    assert section.h == 3.0

def test_rectangular_section_b_setter_accepts_number():
    section = RectangularSection(idt=1, b=2.0, h=3.0)
    section.b = 4.5
    assert section.b == 4.5

def test_rectangular_section_b_setter_rejects_non_number():
    with pytest.raises(TypeError, match="b must be a number"):
        RectangularSection(idt=1, b="x", h=3.0)

def test_rectangular_section_b_setter_rejects_non_positive():
    with pytest.raises(ValueError, match="b must be a positive number"):
        RectangularSection(idt=1, b=0, h=3.0)

def test_rectangular_section_h_setter_accepts_number():
    section = RectangularSection(idt=1, b=2.0, h=3.0)
    section.h = 4.5
    assert section.h == 4.5

def test_rectangular_section_h_setter_rejects_non_number():
    with pytest.raises(TypeError, match="h must be a number"):
        RectangularSection(idt=1, b=2.0, h="x")

def test_rectangular_section_h_setter_rejects_non_positive():
    with pytest.raises(ValueError, match="h must be a positive number"):
        RectangularSection(idt=1, b=2.0, h=0)

def test_rectangular_section_area():
    section = RectangularSection(idt=1, b=2.0, h=3.0)
    assert section.area == 6.0

def test_rectangular_section_perimeter():
    section = RectangularSection(idt=1, b=2.0, h=3.0)
    assert section.perimeter == 10.0

def test_rectangular_section_ix():
    section = RectangularSection(idt=1, b=2.0, h=3.0)
    assert section.ix == pytest.approx(4.5)

def test_rectangular_section_iy():
    section = RectangularSection(idt=1, b=2.0, h=3.0)
    assert section.iy == pytest.approx(2.0)

def test_rectangular_section_it():
    s = RectangularSection(idt=1, b=20, h=60)
    assert s.it == pytest.approx(125308.8454558271)
