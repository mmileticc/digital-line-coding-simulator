from copy import copy

import numpy as np
import matplotlib.pyplot as plt


def generate_binary_sequence(N, P0, P1):
    sequence = np.random.choice([0, 1], size=N, p=[P0, P1])
    return sequence


def unipolar_NRZ(sequence, U):
    return U * sequence


def unipolar_RZ(sequence, U):
    rz_sequence = np.repeat(sequence, 2)
    for i in range(len(rz_sequence)):
        if i % 2 == 1:
            rz_sequence[i] = 0
    return U * rz_sequence


def polar_NRZ(sequence, U):
    return 2 * U * sequence - U


def polar_RZ(sequence, U):
    rz_sequence = np.repeat(sequence, 2)
    rz_sequence = rz_sequence * 2 * U - U
    for i in range(len(rz_sequence)):
        if i % 2 == 1:
            rz_sequence[i] = 0

    return rz_sequence


def manchester(sequence, U):
    manchester_sequence = np.repeat(sequence, 2)
    for i in range(0, len(manchester_sequence), 2):
        if manchester_sequence[i] == 1:
            manchester_sequence[i + 1] = 0
        else:
            manchester_sequence[i + 1] = 1

    return 2 * U * manchester_sequence - U


def AMI(sequence, U):
    ami_sequence = copy(sequence)
    uu = U
    for i in range(len(sequence)):
        if ami_sequence[i] == 1:
            ami_sequence[i] = ami_sequence[i] * uu
            uu *= -1
    return ami_sequence


def m_ary(sequence, U, m, levels):
    b = int(np.log2(m))
    m_ary_sequence = np.zeros(int(len(sequence) / b))

    for i in range(0, len(sequence), b):
        x = 0
        for j in range(b):
            x *= 2
            x += sequence[i + j]
        m_ary_sequence[int(i / b)] = levels[x]
    return m_ary_sequence


def statistika(signal):
    srVrednostSig = np.mean(signal)
    srKvVrSig = (np.mean(np.square(signal)))
    return srVrednostSig, srKvVrSig


def plot_signal(signal, title, sir, U, n=1):
    # plt.bar(range(len(signal)), signal, width=sir, align='edge', color='blue', alpha=0.7)

    sigZaIscrt = np.repeat(signal, n)
    plt.plot(sigZaIscrt)
    # Podesavanje opsega x i y osa
    plt.xlim(0, len(sigZaIscrt) * sir)
    plt.ylim(-U - 0.05, U + 0.05)
    if n < 10:
        n *= 2

    plt.xticks(list(range(0, len(sigZaIscrt) + 1, n)))
    plt.yticks(list(range(-int(U), int(U) + 1, 1)))

    plt.title(title)
    plt.xlabel('Odbirci')
    plt.ylabel('Amplituda')
    plt.grid(True)
    plt.show()


def main():
    N = 300
    n = 10
    U = 1

    levels_4 = [-3 * U, -U, U, 3 * U]
    levels_8 = [-7 * U, -5 * U, -3 * U, -U, U, 3 * U, 5 * U, 7 * U]

    sirinaPodeoka = 1
    titles = ["Unipolarni NRZ", "Unipolarni RZ", "Polarni NRZ", "Polarni RZ",
              "Manchester", "AMI", "M-arni (4 nivoa)", "M-arni (8 nivoa)"]

    ########################### 1. deo pod a) ######################################
    P0 = 0.5
    P1 = 0.5
    # generisanje binarne sekvence
    binary_sequence = generate_binary_sequence(N, P0, P1)

    # generisanje signala svake vrste
    unipolar_NRZ_signal = unipolar_NRZ(binary_sequence, U)
    unipolar_RZ_signal = unipolar_RZ(binary_sequence, U)
    polar_NRZ_signal = polar_NRZ(binary_sequence, U)
    polar_RZ_signal = polar_RZ(binary_sequence, U)
    manchester_signal = manchester(binary_sequence, U)
    AMI_signal = AMI(binary_sequence, U)
    m_ary_4_signal = m_ary(binary_sequence, 4, 4, levels_4)
    m_ary_8_signal = m_ary(binary_sequence, 8, 8, levels_8)

    # racunanje statistika za svaki signal
    signals = [unipolar_NRZ_signal, unipolar_RZ_signal, polar_NRZ_signal, polar_RZ_signal,
               manchester_signal, AMI_signal, m_ary_4_signal, m_ary_8_signal]
    print("-----------------------------------------")
    print("a) Srednje vrednosti kada je P(0) == P(1)")
    print("-----------------------------------------")
    for i in range(len(signals)):
        mean, rms = statistika(signals[i])
        print(f"{titles[i]}:\n\t Srednja vrednost signala = {mean:.3f},\n\t"
              f"Srednja kvadratna vrednost signala= {rms:.3f}")

    ########################### 1. deo pod b) ######################################
    P0 = 0.4
    P1 = 0.6
    # generisanje binarne sekvence
    binary_sequence = generate_binary_sequence(N, P0, P1)

    # generisanje signala svake vrste
    unipolar_NRZ_signal = unipolar_NRZ(binary_sequence, U)
    unipolar_RZ_signal = unipolar_RZ(binary_sequence, U)
    polar_NRZ_signal = polar_NRZ(binary_sequence, U)
    polar_RZ_signal = polar_RZ(binary_sequence, U)
    manchester_signal = manchester(binary_sequence, U)
    AMI_signal = AMI(binary_sequence, U)
    m_ary_4_signal = m_ary(binary_sequence, 4, 4, levels_4)
    m_ary_8_signal = m_ary(binary_sequence, 8, 8, levels_8)

    # racunanje statistika za svaki signal
    signals = [unipolar_NRZ_signal, unipolar_RZ_signal, polar_NRZ_signal, polar_RZ_signal,
               manchester_signal, AMI_signal, m_ary_4_signal, m_ary_8_signal]

    print("----------------------------------------------------")

    print("b) Srednje vrednosti kada je P(1) = 0.6 i P(0) = 0.4")
    print("----------------------------------------------------")
    for i in range(len(signals)):
        mean, rms = statistika(signals[i])
        print(f"{titles[i]}:\n\t Srednja vrednost signala = {mean:.3f},\n\tSrednja kvadratna vrednost signala= {rms:.3f}")

    ########################### 2. deo ######################################

    binary_sequence = generate_binary_sequence(12, P0, P1)
    print("\nDvanaestobitna sekvenca: ", *binary_sequence)

    unipolar_NRZ_signal = unipolar_NRZ(binary_sequence, U)
    unipolar_RZ_signal = unipolar_RZ(binary_sequence, U)
    polar_NRZ_signal = polar_NRZ(binary_sequence, U)
    polar_RZ_signal = polar_RZ(binary_sequence, U)
    manchester_signal = manchester(binary_sequence, U)
    AMI_signal = AMI(binary_sequence, U)
    m_ary_4_signal = m_ary(binary_sequence, 4, 4, levels_4)
    m_ary_8_signal = m_ary(binary_sequence, 8, 8, levels_8)
    signals = [unipolar_NRZ_signal, unipolar_RZ_signal, polar_NRZ_signal, polar_RZ_signal,
               manchester_signal, AMI_signal, m_ary_4_signal, m_ary_8_signal]
    # Prikazivanje oblika signala
    for i in range(len(signals)):
        nn = n
        if i == 1 or i == 3 or i == 4:
            nn /= 2
            nn = int(nn)
        if i == 6:
            plot_signal(signals[i], titles[i], sirinaPodeoka, 4.1, nn)
        elif i == 7:
            plot_signal(signals[i], titles[i], sirinaPodeoka, 8.1, nn)
        else:
            plot_signal(signals[i], titles[i], sirinaPodeoka, U, nn)


if __name__ == '__main__':
    main()
