import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

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

# Параметры для модели распространения COST231
h_BS = 30  # м, высота антенны БС
h_UE = 1.5  # м, высота антенны UE
area_type = 'U'  # тип местности: 'U' - город

# Площади покрытия
Area_open_territory_km2 = 100  # км², площадь открытой территории (макросоты)
Area_indoor_km2 = 4  # км², площадь помещений (микросоты)

# 2. РАСЧЕТ RxSens (ЧУВСТВИТЕЛЬНОСТИ ПРИЕМНИКА)

def calculate_rx_sens(noise_figure_dB, bandwidth_Hz, required_sinr_dB):
    thermal_noise_dBm = -174 + 10 * np.log10(bandwidth_Hz)  # дБм
    rx_sens_dBm = noise_figure_dB + thermal_noise_dBm + required_sinr_dB  # дБм
    return rx_sens_dBm, thermal_noise_dBm

RxSens_BS_dBm, thermal_UL = calculate_rx_sens(NF_BS, BW_UL_Hz, SINR_UL)  # дБм

RxSens_UE_dBm, thermal_DL = calculate_rx_sens(NF_UE, BW_DL_Hz, SINR_DL)  # дБм

print("РАСЧЕТ ЧУВСТВИТЕЛЬНОСТИ ПРИЕМНИКА")

print(f"Тепловой шум (UL, 10 МГц): {thermal_UL:.2f} дБм")
print(f"Тепловой шум (DL, 20 МГц): {thermal_DL:.2f} дБм")
print(f"Чувствительность БС (RxSens_UL): {RxSens_BS_dBm:.2f} дБм")
print(f"Чувствительность UE (RxSens_DL): {RxSens_UE_dBm:.2f} дБм")

# 3. РАСЧЕТ MAPL (МАКСИМАЛЬНО ДОПУСТИМЫХ ПОТЕРЬ)


MAPL_DL = (TxPowerBS - FeederLoss + AntGainBS + MIMOGain -
           IM - PenetrationM) - RxSens_UE_dBm  # дБ

MAPL_UL = (TxPowerUE - FeederLoss + AntGainBS + MIMOGain -
           IM - PenetrationM) - RxSens_BS_dBm  # дБ


print("РАСЧЕТ МАКСИМАЛЬНО ДОПУСТИМЫХ ПОТЕРЬ (MAPL)")

print(f"MAPL для нисходящего канала (DL): {MAPL_DL:.2f} дБ")
print(f"MAPL для восходящего канала (UL): {MAPL_UL:.2f} дБ")

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

d_km_array = np.linspace(0.01, 10, 500)
d_m_array = d_km_array * 1000

pl_fspl = path_loss_fspl(d_km_array, f_GHz)
pl_umi = path_loss_umi_nlos(d_m_array, f_GHz)
pl_cost231 = np.array([path_loss_cost231_hata(d, f_MHz, h_BS, h_UE, area_type) for d in d_km_array])  # дБ

print("МОДЕЛИ РАСПРОСТРАНЕНИЯ СИГНАЛА РАССЧИТАНЫ")

# 5. ПОСТРОЕНИЕ ГРАФИКОВ ЗАВИСИМОСТИ PL(d)

plt.figure(figsize=(12, 8))
plt.plot(d_km_array, pl_fspl, 'g--', label='Free Space (FSPL)', linewidth=1.5)
plt.plot(d_km_array, pl_umi, 'b-', label=f'UMiNLOS (f={f_GHz}ГГц)', linewidth=2)
plt.plot(d_km_array, pl_cost231, 'r-', label=f'COST 231 Hata (Urban, hBS={h_BS}м)', linewidth=2)

plt.axhline(y=MAPL_DL, color='orange', linestyle=':', linewidth=2.5, label=f'MAPL_DL = {MAPL_DL:.1f} дБ')
plt.axhline(y=MAPL_UL, color='purple', linestyle=':', linewidth=2.5, label=f'MAPL_UL = {MAPL_UL:.1f} дБ')

plt.xlabel('Расстояние, d (км)')
plt.ylabel('Потери сигнала, PL (дБ)')
plt.title('Зависимость потерь сигнала от расстояния для разных моделей')
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.legend()
plt.ylim(60, 180)
plt.xlim(0, 10)
plt.savefig('path_loss_models.png', dpi=300, bbox_inches='tight')
plt.show()

print("\nГрафик зависимостей построен и сохранен как 'path_loss_models.png'")

# 6. ОПРЕДЕЛЕНИЕ РАДИУСА СОТЫ

def find_radius(map_loss, distance_array, loss_array):
    interp_func = interp1d(loss_array, distance_array, bounds_error=False, fill_value="extrapolate")
    radius = interp_func(map_loss)
    return radius

radius_DL_cost = find_radius(MAPL_DL, d_km_array, pl_cost231)
radius_UL_cost = find_radius(MAPL_UL, d_km_array, pl_cost231)

radius_DL_umi = find_radius(MAPL_DL, d_m_array, pl_umi)
radius_UL_umi = find_radius(MAPL_UL, d_m_array, pl_umi)

cell_radius_cost = min(radius_DL_cost, radius_UL_cost)
cell_radius_umi = min(radius_DL_umi, radius_UL_umi) / 1000

print("\n" + "=" * 60)
print("ОПРЕДЕЛЕНИЕ РАДИУСА СОТЫ")
print("=" * 60)
print("Для открытой территории (COST231 Hata):")
print(f"  Радиус покрытия DL: {radius_DL_cost:.3f} км")
print(f"  Радиус покрытия UL: {radius_UL_cost:.3f} км")
print(f"  Финальный радиус соты: {cell_radius_cost:.3f} км")

print("\nДля помещений (UMiNLOS):")
print(f"  Радиус покрытия DL: {radius_DL_umi:.1f} м")
print(f"  Радиус покрытия UL: {radius_UL_umi:.1f} м")
print(f"  Финальный радиус соты: {cell_radius_umi:.3f} км")

# 7. РАСЧЕТ ПЛОЩАДИ ПОКРЫТИЯ ОДНОЙ БС

area_per_site_cost = 1.95 * (cell_radius_cost ** 2)
area_per_site_umi = 1.95 * (cell_radius_umi ** 2)


print("РАСЧЕТ ПЛОЩАДИ ПОКРЫТИЯ ОДНОЙ БС")

print(f"Площадь покрытия одной БС (COST231, открытая территория): {area_per_site_cost:.6f} км²")
print(f"Площадь покрытия одной БС (UMiNLOS, помещения): {area_per_site_umi:.6f} км²")

# 8. ВЫЧИСЛЕНИЕ КОЛИЧЕСТВА БАЗОВЫХ СТАНЦИЙ

num_sites_open = np.ceil(Area_open_territory_km2 / area_per_site_cost)  # шт.

num_sites_indoor = np.ceil(Area_indoor_km2 / area_per_site_umi)  # шт.

print("ВЫЧИСЛЕНИЕ КОЛИЧЕСТВА БАЗОВЫХ СТАНЦИЙ")

print(f"Площадь открытой территории: {Area_open_territory_km2} км²")
print(f"Площадь помещений: {Area_indoor_km2} км²")
print(f"\nКоличество БС для открытой территории (COST231): {int(num_sites_open)} шт.")
print(f"Количество БС для помещений (UMiNLOS): {int(num_sites_indoor)} шт.")
print(f"Всего БС для двух территорий: {int(num_sites_open + num_sites_indoor)} шт.")

plt.figure(figsize=(12, 8))
plt.plot(d_km_array, pl_umi, 'b-', label='UMiNLOS модель (помещения)', linewidth=2)
plt.plot(d_km_array, pl_cost231, 'r-', label='COST231 модель (открытая территория)', linewidth=2)

plt.axhline(y=MAPL_UL, color='purple', linestyle='--', linewidth=2, label=f'MAPL_UL = {MAPL_UL:.1f} дБ')
plt.axhline(y=MAPL_DL, color='orange', linestyle='--', linewidth=2, label=f'MAPL_DL = {MAPL_DL:.1f} дБ')

plt.axvline(x=cell_radius_cost, color='red', linestyle=':', linewidth=2,
            label=f'R_открытая = {cell_radius_cost:.2f} км')
plt.axvline(x=cell_radius_umi, color='blue', linestyle=':', linewidth=2,
            label=f'R_помещения = {cell_radius_umi:.3f} км')

plt.xlabel('Расстояние, d (км)')
plt.ylabel('Потери сигнала, PL (дБ)')
plt.title('Определение радиуса соты для разных типов территории')
plt.grid(True, alpha=0.3)
plt.legend()
plt.xlim(0, 6)
plt.ylim(80, 180)
plt.savefig('cell_radius_determination.png', dpi=300, bbox_inches='tight')
plt.show()

print(f"Покрытие открытой территории {Area_open_territory_km2} км² требуется:")
print(f"   {int(num_sites_open)} макросот (модель COST231)")
print(f"   Радиус одной соты: {cell_radius_cost:.2f} км")
print(f"   Площадь покрытия одной БС: {area_per_site_cost:.3f} км²")

print(f"\nПокрытие помещений {Area_indoor_km2} км² требуется:")
print(f"   {int(num_sites_indoor)} микросот (модель UMiNLOS)")
print(f"   Радиус одной соты: {cell_radius_umi:.3f} км ({cell_radius_umi*1000:.1f} м)")
print(f"   Площадь покрытия одной БС: {area_per_site_umi:.4f} км²")

print(f"Общее количество БС для двух территорий: {int(num_sites_open + num_sites_indoor)} шт.")