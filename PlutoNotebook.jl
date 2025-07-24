### A Pluto.jl notebook ###
# v0.19.0

begin
    using PlutoUI
    include(joinpath(@__DIR__, "file.jl"))
end

md"""
# Welcome to the UMERC Pluto Workshop!

This cell includes code from `file.jl` and demonstrates interactivity.
"""

@bind x Slider(1:10, show_value=true)

result = my_function(x)

md"""
**Result from `my_function(x)` in `file.jl`:**

`x = $(x)`, `result = $(result)`
""" 