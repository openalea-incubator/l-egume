from alinea.caribu.sky_tools import GenSky, GetLight, Gensun, GetLightsSun


#soleil
def sunlight():
    soleil=Gensun.Gensun()(1,83,12,43)
    soleil=GetLightsSun.GetLightsSun(soleil)
    soleil=soleil.split(' ')
    soleil=[tuple((float(soleil[0]), tuple((float(soleil[1]), float(soleil[2]), float(soleil[3])))))]

    return soleil



### ciel modele caribu
def determine_sky_Caribu():
    sky = GenSky.GenSky()(1, 'uoc', 3, 4)
    sky_str = GetLight.GetLight(sky)
    #Convert string to list in order to be compatible with CaribuScene input format
    sky_list = []
    for string in sky_str.split('\n'):
        if len(string)!=0:
            string_split = string.split(' ')
            t = tuple((float(string_split[0]), tuple((float(string_split[1]), float(string_split[2]), float(string_split[3])))))
            sky_list.append(t)

    return sky_list

# Calcul poids du ciel Modele de Richard H GRANT
def determine_sky_GRANT(nb_azim, nb_haut, elev_sol, azim_sol):
    tupl_sky = []

    i = 0
    while i < (nb_azim):
        phi = i * (360 / nb_azim)

        j = 0
        while j < nb_haut:
            Elev = j * (90 / nb_haut)
            Zenit_mes = 90 - Elev
            Zenit_sol = 90 - elev_sol

            Psi = np.arccos(
                np.cos(np.radians(Zenit_mes)) * np.cos(np.radians(Zenit_sol)) + np.sin(np.radians(Zenit_mes)) * np.sin(
                    np.radians(Zenit_sol)) * np.cos(np.radians(phi)))

            Energie = 0.0361 * (6.3 + ((1 + np.cos(Psi) * np.cos(Psi)) / (1 - np.cos(Psi)))) * (
                        1 - np.exp(-0.31 / (np.cos(np.radians(Zenit_mes)))))
            dphi = np.radians(360 / nb_azim)
            dteta = np.radians(90 / nb_haut)
            Energie = Energie * np.cos(np.radians(Zenit_mes)) * np.sin(np.radians(Zenit_mes)) * dphi * dteta

            x_cos = np.cos(np.radians(phi)) * np.sin(np.radians(Zenit_mes))
            y_cos = np.sin(np.radians(phi)) * np.sin(np.radians(Zenit_mes))
            z_cos = np.cos(np.radians(Zenit_mes))

            t_sky = tuple((float(Energie), tuple((float(x_cos), float(y_cos), float(-z_cos)))))
            tupl_sky.append(t_sky)

            j += 1

        i += 1
    return tupl_sky



#Calcul poids du ciel Modele du CIE:Spatial distribution of daylight
def determine_sky_CIE(nb_azim, nb_haut, elev_sol, asim_sol, temps_clair):
    tupl_sky = []
    if (temps_clair == True):
        Coeff_a = 4
        Coeff_b = -0.7
        coeff_c = 0
        coeff_d = -1
        coeff_e = 0
    else:
        Coeff_a = -1
        Coeff_b = -0.55
        coeff_c = 2
        coeff_d = -1.5
        coeff_e = 0.15

    i = 0
    while i < (nb_azim):
        asim_mes = i * (360 / nb_azim)

        j = 0
        while j < nb_haut:
            elev_mes = j * (90 / nb_haut)
            Zenit_mes = 90 - elev_mes
            Zenit_sol = 90 - elev_sol

            xi = np.arccos(
                np.cos(np.radians(Zenit_sol)) * np.cos(np.radians(Zenit_mes)) + np.sin(np.radians(Zenit_sol)) * np.sin(
                    np.radians(Zenit_mes)) * np.cos(np.radians(abs(asim_mes - asim_sol))))

            Energie = 1 + coeff_c * (np.exp(coeff_d * xi) - np.exp(coeff_d * np.pi / 2)) + np.exp(1) * np.cos(
                xi) * np.cos(xi)

            x_cos = np.cos(np.radians(asim_mes)) * np.sin(np.radians(Zenit_mes))
            y_cos = np.sin(np.radians(asim_mes)) * np.sin(np.radians(Zenit_mes))
            z_cos = np.cos(np.radians(Zenit_mes))

            t_sky = tuple((float(Energie), tuple((float(x_cos), float(y_cos), float(-z_cos)))))
            tupl_sky.append(t_sky)

            j += 1

        i += 1
    return tupl_sky
