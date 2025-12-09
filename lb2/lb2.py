import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from matplotlib.colors import LinearSegmentedColormap

# 1. ИСХОДНЫЕ ПАРАМЕТРЫ СИСТЕМЫ

# Мощности передачи
TxPowerBS = 46  # дБм, мощность БС
TxPowerUE = 24  # дБм, мощность абонентского устройства
AntGainBS = 21  # дБи, усиление антенны БС
PenetrationM = 15  # дБ, потери на проникновение
IM = 1  # дБ, запас на интерференцию

# Частотные параметры
f_MHz = 1800  # МГц, рабочая частота
f_GHz = f_MHz / 1000  # ГГц, рабочая частота
BW_UL_Hz = 10e6  # Гц, полоса UL канала
BW_DL_Hz = 20e6  # Гц, полоса DL канала

# Шумовые коэффициенты
NF_BS = 2.4  # дБ, коэффициент шума БС
NF_UE = 6  # дБ, коэффициент шума UE

# Требования к качеству связи
SINR_DL = 2  # дБ, минимальный SINR для DL
SINR_UL = 4  # дБ, минимальный SINR для UL

# Технологические параметры
MIMOGain = 3  # дБ, усиление MIMO
FeederLoss = 2.9  # дБ, потери в фидерном тракте (2 + 0.5 + 0.4)

# Параметры для модели распространения
h_BS = 30  # м, высота антенны БС
h_UE = 1.5  # м, высота антенны UE
area_type = 'U'  # тип местности: 'U' - город

# Параметры для Walfish-Ikegami
h_roof = 20  # м, средняя высота крыш зданий
w_street = 20  # м, ширина улиц
phi = 90  # градусов, угол между направлением связи и осью улицы
b = 30  # м, расстояние между зданиями


# 2. РАСЧЕТ RxSens (ЧУВСТВИТЕЛЬНОСТИ ПРИЕМНИКА)

def calculate_rx_sens(noise_figure_dB, bandwidth_Hz, required_sinr_dB):
    thermal_noise_dBm = -174 + 10 * np.log10(bandwidth_Hz)  # дБм
    rx_sens_dBm = noise_figure_dB + thermal_noise_dBm + required_sinr_dB  # дБм
    return rx_sens_dBm, thermal_noise_dBm


RxSens_BS_dBm, thermal_UL = calculate_rx_sens(NF_BS, BW_UL_Hz, SINR_UL)  # дБм
RxSens_UE_dBm, thermal_DL = calculate_rx_sens(NF_UE, BW_DL_Hz, SINR_DL)  # дБм

# 3. РАСЧЕТ MAPL (МАКСИМАЛЬНО ДОПУСТИМЫХ ПОТЕРЬ)

MAPL_DL = (TxPowerBS - FeederLoss + AntGainBS + MIMOGain -
           IM - PenetrationM) - RxSens_UE_dBm  # дБ

MAPL_UL = (TxPowerUE - FeederLoss + AntGainBS + MIMOGain -
           IM - PenetrationM) - RxSens_BS_dBm  # дБ


# 4. МОДЕЛИ РАСПРОСТРАНЕНИЯ СИГНАЛА

def path_loss_fspl(d_km, f_GHz):
    d_m = d_km * 1000  # м
    pl_dB = 20 * np.log10((4 * np.pi * d_m * f_GHz * 1e9) / 3e8)  # дБ
    return pl_dB


def path_loss_umi_nlos(d_m, f_GHz):
    pl_dB = 26 * np.log10(f_GHz) + 22.7 + 36.7 * np.log10(d_m)  # дБ
    return pl_dB


def path_loss_cost231_hata(d_km, f_MHz, h_BS, h_UE, area_type):
    if 150 <= f_MHz <= 1500:
        A, B = 69.55, 26.16
    elif 1500 < f_MHz <= 2000:
        A, B = 46.3, 33.9

    if area_type in ['DU', 'U']:
        a = 3.2 * (np.log10(11.75 * h_UE)) ** 2 - 4.97
    elif area_type in ['SU', 'RURAL', 'ROAD']:
        a = (1.1 * np.log10(f_MHz) - 0.7) * h_UE - (1.56 * np.log10(f_MHz) - 0.8)

    if area_type == 'DU':
        L_clutter = 3  # дБ
    elif area_type == 'U':
        L_clutter = 0  # дБ
    elif area_type == 'SU':
        L_clutter = -(2 * (np.log10(f_MHz / 28)) ** 2 + 5.4)  # дБ
    elif area_type == 'RURAL':
        L_clutter = -(4.78 * (np.log10(f_MHz)) ** 2 - 18.33 * np.log10(f_MHz) + 40.94)  # дБ
    elif area_type == 'ROAD':
        L_clutter = -(4.78 * (np.log10(f_MHz)) ** 2 - 18.33 * np.log10(f_MHz) + 35.94)  # дБ
    else:
        L_clutter = 0  # дБ

    if d_km >= 1:
        s = 44.9 - 6.55 * np.log10(h_BS)
    else:
        s = (47.88 + 13.9 * np.log10(f_MHz) - 13.9 * np.log10(h_BS)) * (1 / np.log10(50))

    pl_dB = (A + B * np.log10(f_MHz) - 13.82 * np.log10(h_BS) - a +
             s * np.log10(d_km) + L_clutter)  # дБ
    return pl_dB


def path_loss_walfish_ikegami(d_km, f_MHz, h_BS, h_UE, h_roof, w_street=20, phi=90, b=30):

    d_m = d_km * 1000
    f = f_MHz

    L0 = 32.44 + 20 * np.log10(f) + 20 * np.log10(d_m)

    if d_m <= 100 and h_BS > h_roof + 5:
        L_los = 42.6 + 26 * np.log10(d_m) + 20 * np.log10(f/1000)
        return L_los

    if 0 <= phi < 35:
        Lori = 0
    elif 35 <= phi < 55:
        Lori = 2.5 + 0.075 * (phi - 35)
    else:
        Lori = 4.0 - 0.114 * (phi - 55)

    Lrts = -16.9 - 10 * np.log10(w_street) + 10 * np.log10(f) + \
           20 * np.log10(h_roof - h_UE) + Lori

    if h_BS > h_roof:
        Lbsh = 0
        ka = 54
        kd = 18
    else:
        delta_h = h_roof - h_BS
        if delta_h > 0:
            Lbsh = -18 * np.log10(1 + (h_BS - h_roof))
            if d_m >= 500:
                ka = 54 - 0.8 * (h_BS - h_roof)
            else:
                ka = 54 - 1.6 * (h_BS - h_roof) * d_m / 1000
            kd = 18 - 15 * (h_BS - h_roof) / h_roof
        else:
            Lbsh = 0
            ka = 54
            kd = 18

    if f >= 1500:
        kf = -4 + 0.7 * (f / 925 - 1)
    else:
        kf = -4 + 1.5 * (f / 925 - 1)

    Lmsd = Lbsh + ka + kd * np.log10(d_m) + kf * np.log10(f) - 9 * np.log10(b)

    L_total = L0 + Lrts + Lmsd

    return L_total


d_km_array = np.linspace(0.01, 10, 500)
d_m_array = d_km_array * 1000

pl_fspl = path_loss_fspl(d_km_array, f_GHz)
pl_umi = path_loss_umi_nlos(d_m_array, f_GHz)
pl_cost231 = np.array([path_loss_cost231_hata(d, f_MHz, h_BS, h_UE, area_type) for d in d_km_array])
pl_walfish = np.array([path_loss_walfish_ikegami(d, f_MHz, h_BS, h_UE, h_roof, w_street, phi, b) for d in d_km_array])

# 5. ПОСТРОЕНИЕ ГРАФИКА

plt.figure(1)
plt.plot(d_km_array, pl_fspl, 'g--', label='Free Space (FSPL)', linewidth=1.5)
plt.plot(d_km_array, pl_umi, 'b-', label=f'UMiNLOS (f={f_GHz}ГГц)', linewidth=2)
plt.plot(d_km_array, pl_cost231, 'r-', label=f'COST 231 Hata (Urban, hBS={h_BS}м)', linewidth=2)
plt.plot(d_km_array, pl_walfish, 'm-', label=f'Walfish-Ikegami (h_roof={h_roof}м)', linewidth=2)

plt.axhline(y=MAPL_DL, color='orange', linestyle=':', linewidth=2.5, label=f'MAPL_DL = {MAPL_DL:.1f} дБ')
plt.axhline(y=MAPL_UL, color='purple', linestyle=':', linewidth=2.5, label=f'MAPL_UL = {MAPL_UL:.1f} дБ')

plt.xlabel('Расстояние, d (км)')
plt.ylabel('Потери сигнала, PL (дБ)')
plt.title('Зависимость потерь сигнала от расстояния для разных моделей (с Walfish-Ikegami)')
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.legend()
plt.ylim(60, 180)
plt.xlim(0, 10)
plt.savefig('path_loss_models_with_walfish.png', dpi=300, bbox_inches='tight')
plt.show()

print("График с моделью Walfish-Ikegami сохранен как 'path_loss_models_with_walfish.png'")

# 6. ТЕПЛОВАЯ КАРТА РАСПРОСТРАНЕНИЯ СИГНАЛА ДЛЯ WALFISH-IKEGAMI

max_radius_km = 6.0

x = np.linspace(-max_radius_km, max_radius_km, 200)
y = np.linspace(-max_radius_km, max_radius_km, 200)
X, Y = np.meshgrid(x, y)

R = np.sqrt(X ** 2 + Y ** 2)

EIRP = TxPowerBS - FeederLoss + AntGainBS + MIMOGain



signal_strength = np.zeros_like(R)
for i in range(len(x)):
    for j in range(len(y)):
        distance_km = R[j, i]
        if distance_km > 0:
            pl = path_loss_walfish_ikegami(distance_km, f_MHz, h_BS, h_UE, h_roof, w_street, phi, b)
            rx_power = EIRP - pl
            signal_strength[j, i] = rx_power
        else:
            signal_strength[j, i] = EIRP


# Находим радиус покрытия для Walfish-Ikegami
def find_radius(map_loss, distance_array, loss_array):
    interp_func = interp1d(loss_array, distance_array, bounds_error=False, fill_value="extrapolate")
    radius = interp_func(map_loss)
    return radius


radius_DL_walfish = find_radius(MAPL_DL, d_km_array, pl_walfish)
radius_UL_walfish = find_radius(MAPL_UL, d_km_array, pl_walfish)
cell_radius_walfish = min(radius_DL_walfish, radius_UL_walfish)

print(f"\nРадиус соты для Walfish-Ikegami: {cell_radius_walfish:.3f} км")


plt.figure(2)

heatmap = plt.contourf(X, Y, signal_strength, levels=30, cmap='RdYlBu_r')
plt.colorbar(heatmap, label='Мощность сигнала, дБм')
plt.scatter(0, 0, color='red', s=200, label='Базовая станция', edgecolors='black')
plt.contour(X, Y, R, levels=[cell_radius_walfish], colors='green', linewidths=2, linestyles='--')
plt.xlabel('Расстояние по оси X, км')
plt.ylabel('Расстояние по оси Y, км')
plt.title('Тепловая карта распространения сигнала (Walfish-Ikegami)')
plt.legend()
plt.grid(True, alpha=0.3)
plt.axis('equal')

radius_DL_cost = find_radius(MAPL_DL, d_km_array, pl_cost231)
radius_UL_cost = find_radius(MAPL_UL, d_km_array, pl_cost231)
cell_radius_cost = min(radius_DL_cost, radius_UL_cost)

plt.figure(3)
plt.axvline(x=cell_radius_cost, color='red', linestyle=':', linewidth=2,
            label=f'R_COST231 = {cell_radius_cost:.2f} км')
plt.axvline(x=cell_radius_walfish, color='magenta', linestyle=':', linewidth=2,
            label=f'R_Walfish = {cell_radius_walfish:.2f} км')

plt.xlabel('Расстояние, км')
plt.ylabel('Потери сигнала, дБ')
plt.title('Сравнение COST231 и Walfish-Ikegami')
plt.legend()
plt.grid(True, alpha=0.3)
plt.xlim(0, max(cell_radius_cost, cell_radius_walfish) * 1.5)

plt.tight_layout()
plt.show()

