S = S4.NewSimulation()
pcall(loadstring(S4.arg))

S:SetLattice({period, 0}, {0, 0})

S:SetNumG(harmonics)

S:AddMaterial(
    'Cover',
    {
        (cover_n * cover_n) - (cover_k * cover_k),
        (2 * cover_n * cover_k)
    }
) 
S:AddMaterial(
    'Grating',
    {
        (material_n * material_n) - (material_k * material_k),
        (2 * material_n * material_k)
    }
)
S:AddMaterial(
    'Waveguide',
    {
        (material_n * material_n) - (material_k * material_k),
        (2 * material_n * material_k)
    }
)
S:AddMaterial(
    'Substrate',
    {
        (substrate_n * substrate_n) - (substrate_k * substrate_k),
        (2 * substrate_n * substrate_k)
    }
)

waveguide_thickness = film_thickness - grating_thickness + 0.01
S:AddLayer('cover', 0, 'Cover')
S:AddLayer('grating', grating_thickness, 'Grating')
S:AddLayer('waveguide', waveguide_thickness, 'Waveguide')
S:AddLayer('substrate', 0, 'Substrate')

S:SetLayerPatternRectangle(
	'grating',
	'Cover',
	{0, 0},
	0,
	{0.5 * (1 - (fill_factor)) * period, 0}
)

S:SetExcitationPlanewave(
	{0, 0},
        {TE, 0},
        {TM, 0})

for lambda = wavelength_initial, wavelength_final, wavelength_step do

	freq = 1/lambda
	S:SetFrequency(freq)

	transmission = S:GetPowerFlux('substrate')
	inc, reflection = S:GetPowerFlux('cover', 10)

	print(lambda, transmission, -reflection)
end
