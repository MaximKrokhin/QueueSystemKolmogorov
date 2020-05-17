import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)
import math
import time
start = time.time()
def rungeKutta(f, t0, p0, tEnd, tau):
         def increment(f, t, p, tau):# поиск приближённого решения методом Рунге—Кутта—Фельберга.
                  k1 = tau * f(t, p)
                  k2 = tau * f(t + (1/4) * tau, p + (1/4) * k1)
                  k3 = tau *f(t + (3/8) * tau, p + (3/32) * k1 + (9/32) * k2)
                  k4 = tau * f(t + (12/13) * tau, p + (1932/2197) * k1 - (7200/2197) * k2 + (7296/2197) * k3)
                  k5 = tau*f(t + tau, p + (439/216) * k1 - 8 * k2 + (3680/513) * k3 - (845/4104) * k4)
                  k6 = tau*f(t + (1/2) * tau, p - (8/27) * k1 + 2 * k2 - (3544/2565) * k3 + (1859/4104) * k4 - (11/40) * k5)
                  return (16/135) * k1 + (6656/12825) * k3 + (28561/56430) * k4 - (9/50) * k5 + (2/55) * k6
         t = []#подготовка пустого списка t
         p = []#подготовка пустого списка p
         t.append(t0)#внесение в список t начального значения to
         p.append(p0)#внесение в список p начального значения po
         while t0 < tEnd:#внесение результатов расчёта в массивы t,p
                  tau = min(tau, tEnd - t0)#определение минимального шага tau
                  p0 = p0 + increment(f, t0, p0, tau) # расчёт значения в точке t0,p0 для задачи Коши
                  t0 = t0 + tau # приращение времени
                  t.append(t0) # заполнение массива t
                  p.append(p0) # заполнение массива p
         return np.array(t),np.array(p)

def f(t, p): #система уравнений
   func = np.zeros([n+m+1])
   summ= 1-sum(p)+p[n+m]
   func[0] = -1 * lam * p[0] + n * mu * p[1]
   for k in range(1, n+m-1):
      func[k] = lam * p[k-1] - (lam + n * mu) * p[k] + n * mu * p[k+1]
   func[n + m - 1] = lam * p[n + m - 2] - (lam + n * mu) * p[n + m -1] + n * mu * summ
   func[n + m] = lam * p[n + m - 1] - n * mu * summ
   return func

n = 5 #начальные параметры - число каналов
m = 6 #места в очереди
t0 = 0.0 #начальный момент времени
tEnd = 50.0 #конечный момент времени
p0 = np.zeros([n+m+1]) #начальные вероятности
p0[0] = 1.0 #начальное значение вероятности состояние 0 заявок
lam = 0.5 #интенсивность поступления заявок
mu = 0.15 #интенсивность обслуживания заявок
tau = 0.01 #величина шага
hi = lam/(n*mu) #коэф-т загруженности

A = lam*(1-hi**(n+m))/(1-hi**(n+m+1))
print("Абсолютная пропускная способность: ",A)
q = (1-hi**(n+m))/(1-hi**(n+m+1))
print("Относительная пропускная способность: ", q)
p_rej = (1-hi)*(hi**(n+m))/(1 - hi**(m+n+1))
print("Вероятность отказа: ", p_rej)
p_hold = 1 -((1-hi**m)*(hi**(n+1))/(1-hi**(n+m+1)))
print("Вероятность отсутствия очереди: ", p_hold)
l = (1-(n+m+1)*(hi**(n+m))+(n+m)*(hi**(n+m+1)))*hi/((1-hi)*(1-hi**(n+m+1)))
print("Среднее число заявок в системе: ", l)
t_serv_avg = ((1-(n+1)*(hi**n)+n*hi**(n+1))/(n*mu*(1-hi)*(1-hi**(n+m+1))))+(((hi**n)*(1-hi**m))/(mu*(1-hi**(n+m+1))))
print("Среднее время обслуживания заявки в системе: ", t_serv_avg)
r = (1-(m+1)*(hi**m)+m*(hi**(m+1)))*hi/((1-hi)*(1-hi**(n+m+1)))
print("Среднее число заявок в очереди: ", r)
t_queue_avg = (1 - (m+1)*(hi**m)+m*(hi**(m+1)))/(n*mu*(1-hi)*(1-hi**(n+m+1)))
print("Среднее время заявки в очереди", t_queue_avg)
t_smo_avg = (1 - (n+m+1)*(hi**(n+m))+(n+m)*(hi**(n+m+1)))/(n*mu*(1-hi)*(1-hi**(n+m+1)))
print("Среднее время нахождения заявки в системе", t_smo_avg)
t, p = rungeKutta(f, t0, p0, tEnd, tau)

#интенсивность потока
pow_hi = np.array([hi**i for i in range(0,n+m+1)])

#вектор предельных вероятностей
p_lim=list()
sum_hi = np.sum(pow_hi)
p_lim.append(1/sum_hi)
for i in range (1,n+m+1):
    p_lim.append(pow_hi[i]/sum_hi)
print(p_lim, sum(p_lim))

stop = time.time()

print ("Время на модельную задачу: %f"%(stop-start))
fig, ax = plt.subplots(figsize=(10, 8))
plt.plot(t, p)  # графики
plt.xlabel('Time')
plt.ylabel('Probability')
ax.set_xlim(0, 50)
ax.set_ylim(0, 1)
ax.xaxis.set_major_locator(MultipleLocator(2))
ax.yaxis.set_major_locator(MultipleLocator(0.05))
ax.xaxis.set_minor_locator(AutoMinorLocator(1))
ax.yaxis.set_minor_locator(AutoMinorLocator(0.01))
ax.grid(which='major', color='#CCCCCC', linestyle='--')
ax.grid(which='minor', color='#CCCCCC', linestyle=':')
plt.grid(True)
plt.show()


