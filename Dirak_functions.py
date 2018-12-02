# ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ

import math
import scipy
from scipy import integrate
import numpy

# Функция Ферми-Дирака I_1/2
def integral_1_2(x):
    sqx = x ** (1 / 2)
    xsqx = x * sqx
    # Апроксимация функции дирака трёхчленной формулой при k=1/2
    if (x >= 176868.709):
        return xsqx / 1.5 + 0.82246703342411321823620758332301 / sqx + 0.71027462212293443818237742585514 / (x * xsqx);
    elif (x < -17.9):
        return 0.88622692545275801364908374167057 * math.exp(x)
    else:
        if (x > 30.0):
            i0 = x
        else:
            i0 = math.log(1.0 + math.exp(x), math.e)

        return (0.8862269254528) * i0 * ((1.0 + 1.178 * i0 + (0.1812102675) * (i0 ** 3.0)) ** (1.0 / 6.0))

    # Функция Ферми-Дирака I_3/2


def integral_3_2_1(x):
    return 3 / 10 * integral_1_2(x) * (125 + 60 * integral_1_2(x) + 18 * (integral_1_2(x)) ** 2) ** (1 / 3)



def integral_3_2(x):

    if (x >= 10.3392439415065):
        cf0 = 0.4
        cf1 = 2.467401100272339654708622749969

        cf2 = -0.71027462212293443818237742585514

        cf3 = -2.7718624442740362302768603751139

        cf4 = -44.13000875394151241108382894273

        cf5 = -15821.16373906104901524376218556

        cf6 = -100930.88369938006505373290376014

        sqx = x**(1/2)
        xsqx = x * sqx
        return cf0 * x * xsqx + cf1 * sqx + cf2 / xsqx + cf3 * x ** (-3.5) + cf4 * x ** (-5.5) + cf5 * x ** (-7.5) + cf6 * x ** (-9.5)


    if(x >= 4.0):

        kFac = 3.23171853624284
        bFac = 7.16962197075325
        cf0 = 7.19740367036291E1
        cf1 = 6.47906294576754E1
        cf2 = 1.0339089561051E1
        cf3 = 4.11667482183952E-1
        cf4 = -1.37618462794111E-2
        cf5 = 1.15494755799039E-3
        cf6 = -1.2412400150752E-4
        cf7 = 1.28217367911798E-5
        return T7_Cheb(x, kFac, bFac, cf0, cf1, cf2, cf3, cf4, cf5, cf6, cf7)
    if(x >= 1.0):

        kFac = 1.52938673731248
        bFac = 2.5
        cf0 = 8.95753044162565E0
        cf1 = 7.56437857186623E0
        cf2 = 1.28658303215612E0
        cf3 = 7.69768485869776E-2
        cf4 = -2.11582983450359E-3
        cf5 = -4.39693320860579E-5
        cf6 = 2.74602973946436E-5
        cf7 = -3.60319354575367E-6
        return T7_Cheb(x, kFac, bFac, cf0, cf1, cf2, cf3, cf4, cf5, cf6, cf7)
    if(x >= -0.994315794013506):

        kFac = 1.00979557910416
        bFac = 9.79557910415929E-3
        cf0 = 1.37154508993631E0
        cf1 = 1.09969157706526E0
        cf2 = 2.09559595052497E-1
        cf3 = 2.14341276992499E-2
        cf4 = 8.1511711708937E-4
        cf5 = -6.0503775127313E-5
        cf6 = -5.3080886762355E-6
        cf7 = 5.55544635783695E-7
        return T7_Cheb(x, kFac, bFac, cf0, cf1, cf2, cf3, cf4, cf5, cf6, cf7)
    else:
        expx = math.exp(x)
        cf0 = 1.3293403881791370204736256125059

        cf1 = 0.17677669529663688110021109052621
        cf2 = 0.06415002990995841827879430894466
        cf3 = 0.03125
        cf4 = 0.01788854381999831757127338934985
        cf5 = 0.01134023029066286156572816701253
        cf6 = 0.00771356067365769851458197012723
        cf7 = 0.00552427172801990253438159657894
        cf8 = 0.00411522633744855967078189300412
        cf9 = 0.00316227766016837933199889354443
        return cf0 * expx * (
                    1.0 - cf1 * expx + cf2 * math.exp(2 * x) - cf3 * math.exp(3 * x) + cf4 * math.exp(4 * x) - cf5 * math.exp(5 * x) + cf6 * math.exp(6 * x) - cf7 * math.exp(7 * x) + cf8 * math.exp(8 * x) - cf9 * math.exp(9 * x))



def integral_minus_1_2(x):
    if (x > 500000.0):
        return 2 * x ** (1 / 2)
    else:
        if (x < -17.9):
            return 2 * 0.88622692545275801364908374167057 * math.exp(x)
        else:
            if (x > 30.0):
                i0 = x
            else:
                i0 = math.log(1.0 + math.exp(x))

            return 2 * (0.8862269254528) * i0 * ((1.0 + 1.614 * i0 + (0.4844730731296) * (i0 ** 3.0)) ** (-1.0 / 6.0))

def igrek_sht(x):
    dx = max(0.14, 0.001 * abs(x))
    return (-1/2*integral_minus_1_2(x + dx) + 1/2*integral_minus_1_2(x - dx) + 168 * 1/2*integral_minus_1_2(x + dx / 2.0) - 168 * 1/2*integral_minus_1_2(x - dx / 2.0) - 5376 * 1/2*integral_minus_1_2(x + dx / 4.0) + 5376 * 1/2*integral_minus_1_2(x - dx / 4.0) + 32768 * 1/2*integral_minus_1_2(x + dx / 8.0) - 32768 * 1/2*integral_minus_1_2(x - dx / 8.0)) / (5670 * dx)



#Сумма полиномов Чебышева
def T7_Cheb(x,kFac,bFac,cf0,cf1,cf2,cf3,cf4,cf5,cf6,cf7):
    xT=(x-bFac)/kFac
    x2T=xT*xT
    x3T=xT**3
    x4T=xT**4
    x5T=xT**5
    x6T=xT**6
    x7T=xT**7
    return cf0+cf1*xT+cf2*(2.0*x2T-1.0)+cf3*(4.0*x3T-3.0*xT)+cf4*(8.0*x4T-8.0*x2T+1.0)+cf5*(16.0*x5T-20.0*x3T+5.0*xT)+cf6*(32.0*x6T-48.0*x4T+18.0*x2T-1.0)+cf7*(64.0*x7T-112.0*x5T+56.0*x3T-7.0*xT)


#Аппроксимация только обменного интеграла(без членов для ТФП)

def J_Exchange(x):

    if (x >= 15.0):
        bFac = 0.383708729538668
        cf1 = -0.82246703342411321823620758332301

        cf2 = 0.5
        cf3 = 1.691130052673653424243755775845575

        cf4 = 14.18716691908001339367511331778745

        cf5 = 346.0733028286156983792120782023865
        cf6 = 16688.67955429090163804075587296032
        cf7 = 1.336392345120319416764032215294369*10**8
        cf8 = 1.603501910160424576402142699306101*10**8
        bFac = 0.383705176155388
        x2 = x * x
        return  bFac + cf1 * math.log(x) + cf2 * x2 + cf3 / x2 + cf4 * x ** (-4) + cf5 * x ** (-6) + cf6 * x ** (-8)
    if (x >= 7.0):

        kFac = 4.07836463283327
        bFac = 11
        cf0 = 6.31198547629997E1
        cf1 = 4.45291937700006E1
        cf2 = 4.1939724684287E0
        cf3 = -5.37494213520073*10**(-3)
        cf4 = 9.20354291643832*10**(-4)
        cf5 = -1.5956682533913*10**(-4)
        cf6 = 2.53521759603359*10**(-5)
        cf7 = -3.28663978321764*10**(-6)
        return T7_Cheb(x, kFac, bFac, cf0, cf1, cf2, cf3, cf4, cf5, cf6, cf7)
    if(x >= 3.0):

        kFac = 2.03918231641664
        bFac = 5
        cf0 = 1.27326204335739E1
        cf1 = 9.78231278110018E0
        cf2 = 1.08511507707662E0
        cf3 = -3.81689002985508E-3
        cf4 = -2.54316506324681E-4
        cf5 = 1.91871819682721E-4
        cf6 = -4.6674660811763E-5
        cf7 = 7.50800003121444E-6
        return T7_Cheb(x, kFac, bFac, cf0, cf1, cf2, cf3, cf4, cf5, cf6, cf7)
    if(x >= 1.0):

        kFac = 1.01959115820832
        bFac = 2
        cf0 = 2.20669406679466E0
        cf1 = 1.73765882478791E0
        cf2 = 2.46276242685701E-1
        cf3 = 6.95198753002888E-3
        cf4 = -9.09023628906075E-4
        cf5 = 3.99027652310952E-5
        cf6 = 5.61013783821618E-6
        cf7 = -1.10487782678526E-6
        return T7_Cheb(x, kFac, bFac, cf0, cf1, cf2, cf3, cf4, cf5, cf6, cf7)
    if(x >= 0.0):
        kFac = 0.509795579104159
        bFac = 0.5
        cf0 = 4.27168174319067E-1
        cf1 = 2.67313106138003E-1
        cf2 = 3.51723157784946E-2
        cf3 = 2.0163929819443E-3
        cf4 = -4.37501397275442E-6
        cf5 = -5.25135086302045E-6
        cf6 = 5.26130988671991E-8
        cf7 = 2.20201508919354E-8
        return T7_Cheb(x, kFac, bFac, cf0, cf1, cf2, cf3, cf4, cf5, cf6, cf7)
    if(x >= -1.48387892035288):

        kFac = 0.757346684328119
        bFac = -0.742653315671881
        cf0 = 8.41142629137387E-2
        cf1 = 8.82401311543843E-2
        cf2 = 2.36549814019604E-2
        cf3 = 3.73539718344219E-3
        cf4 = 3.24712588175925E-4
        cf5 = 4.40494239728169E-6
        cf6 = -2.30871914634807E-6
        cf7 = 1.97582995423813E-7
        return T7_Cheb(x, kFac, bFac, cf0, cf1, cf2, cf3, cf4, cf5, cf6, cf7)
    else:
        exp2x = math.exp(2 * x)
        cf0 = 0.78539816339744830961566084581988
        cf1 = 0.47140452079103168293389624140323
        cf2 = 0.41367513459481288225457439025098
        cf3 = 0.36329931618554520654648560498039
        cf4 = 0.32247788425329945604964119382178
        cf5 = 0.28947176887871823766305078775578
        cf6 = 0.26245962433780030878692715495100
        cf7 = 0.24002748681137621288729803322376
        cf8 = 0.22113507376025144455650756499489
        cf9 = 0.20502010799708259411082291098565
        cf10 = 0.19111818527221522950993007500165
        return cf0 * exp2x * (0.5 - cf1 * math.exp(x) + cf2 * exp2x - cf3 * math.exp(3 * x) + cf4 * math.exp(4 * x) - cf5 * math.exp(5 * x) + cf6 * math.exp(6 * x) - cf7 * math.exp(7 * x) + cf8 * math.exp(8 * x) - cf9 * math.exp(9 * x) + cf10 * math.exp(10 * x))


def igrek(x):
    return 6*J_Exchange(x)+integral_1_2(x)*1/2*integral_minus_1_2(x)