
function [ p_dB_SPL ] = convert_complex_pressure_to_dB_SPL( complex_pressures )


p_mag = abs( complex_pressures );
    p_rms = p_mag / sqrt(2);
        p_dB_SPL = 20*log10( p_rms / 20e-6 );


