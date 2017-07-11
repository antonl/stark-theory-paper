import pint

q = pint.UnitRegistry()
D = 3.34e-30*q('C m')
A3 = 1.113e-40*q('C m^2 V^-1')
to_wn = lambda x: (x/(q.planck_constant*q.speed_of_light)).to('1/cm')

def convert_temperature(kelvin):
    '''convert temperature to 1/cm (energy) units used in simulations
    
    :param kelvin: temperature in Kelvin
    :return: 
    '''
    return to_wn(q.boltzmann_constant*kelvin*q('K'))

if __name__ == '__main__':
    print("77K: ", convert_temperature(77))
