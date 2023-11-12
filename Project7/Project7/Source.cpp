#include<iostream>
#include<array>
#include<cmath>
#include<iomanip>


const unsigned int N = 5;
std::array <double, 15> h; //массив узлов

template<typename RealType, unsigned int N> //структура для вывода
struct DerivativeCoef
{
	RealType centralCoef;
	std::array<RealType, N> otherCoefs;
};

//Реализация метода Гаусса
template<typename RealType, unsigned int N>
DerivativeCoef<RealType, N> GaussMethod(std::array<std::array<RealType, N + 2>, N + 1>& x)
{
	RealType centralCoefs, dub, dub2, maxElem, sum = 0;
	int maxIndex;
	std::array<RealType, N> otherCoefs;


	for (int(l) = 0; l < N + 1; l++)
	{
		maxIndex = l;
		maxElem = x[l][l]; //изначально считаем, что максимальный элемент верхний левый
		for (int(i) = l + 1; i < N + 1; i++) //ищем максимальный элемент в столбце, чтобы переставить строку наверх
		{
			if (abs(x[i][l]) > abs(maxElem))
			{
				maxIndex = i;
				maxElem = x[i][l]; //должен быть типа RealType
			}
		}
		//меняем строки местами
		for (int(j) = l; j < N + 2; j++)
		{
			dub = x[l][j]; //должен быть типа RealType
			x[l][j] = x[maxIndex][j];
			x[maxIndex][j] = dub;

		}
		//нормируем строки по первому элементу
		for (int(i) = l; i < N + 1; i++)
		{
			dub2 = x[i][l]; //Тип RealType, дублер первого элемента в каждой строке, на который делим всю строку
			for (int(j) = l; j < N + 2; j++)
			{
				if (dub2 != 0)
				{
					x[i][j] = x[i][j] / dub2;
				}

				if (i != l && dub2 != 0)
				{
					x[i][j] = x[i][j] - x[l][j]; //вычитание строк одной из другой
				}
			}

		}


	}

	//нашли последний корень
	otherCoefs[N - 1] = x[N][N + 1];
	//ищем остальные корни
	for (int(i) = N - 1; i >= 0; i--)
	{
		for (int(j) = i + 1; j < N + 1; j++)
		{
			sum = sum + x[i][j] * otherCoefs[j - 1]; //sum типа RealType 
		}
		if (i == 0)
		{
			centralCoefs = x[i][N + 1] - sum;
		}
		else
		{
			otherCoefs[i - 1] = x[i][N + 1] - sum;
			sum = 0;
		}
	}

	DerivativeCoef<RealType, N> S;
	S.centralCoef = centralCoefs;
	S.otherCoefs = otherCoefs;

	return S;
}



template<typename RealType, unsigned int N>
DerivativeCoef<RealType, N> calcDerivativeCoef(const std::array<RealType, N>& points) noexcept //функция, которая считает коэффициенты производные
{

	std::array<std::array<double, N + 2>, N + 1> x; //матрица коэффициентов для слау
	for (int(i) = 0; i < N + 1; i++)
	{
		for (int(j) = 0; j < N + 2; j++)
		{
			if (i == 0) //рассматриваем первую строку матрицы
			{
				if (j != N + 1) //в первой строке все элементы (не включая последний) равны 1
					x[i][j] = 1;
				else x[i][j] = 0; // в первой строке последний элемент равен 0
			}
			if (i != 0) //рассматриваем остальные строки после первой
			{
				if (j == 0) // в первом столбце все элементы равны нулю, так как рассматриваем центральную точку
				{
					x[i][j] = 0;
				}
				if (j != 0 && j != N + 1) //вычисление из разложения в ряд Тейлора всех элементов кроме последнего
				{
					x[i][j] = pow(points[j - 1], i);
				}
				if (j == N + 1 && i == 1) //последний элемент во второй строке равен 1
				{
					x[i][j] = 1;
				}
				if (j == N + 1 && i != 1) //последний элемент в остальных строках равен 0
				{
					x[i][j] = 0;
				}
			}
		}
	}
	//вывод матрицы
	std::cout << "Матрица коэффициентов для СЛАУ:" << std::endl;
	for (int(i) = 0; i < N + 1; i++) 
	{
		for (int(j) = 0; j < N + 2; j++)
		{
			std::cout << x[i][j] << " "; //массив коэффициентов для узлов
		}
		std::cout << std::endl;
	}

	return GaussMethod<RealType, N>(x);

}

////////////////////////////////////////////
//////////////////////////////MAIN/////////

int main()
{
	setlocale(LC_ALL, "Russian");

	std::cout << "N = "<<N << std::endl;

	std::array <double, N> points; //массив узлов
	points = { -2,-1,1,2,3 };

	std::cout << "Узлы:" << std::endl; //вывод массива узлов
	for (int(i) = 0; i < N; i++)
	{
		std::cout << points[i] << " ";
	}
	std::cout << std::endl;
	DerivativeCoef<double, N> coefs; //объект структуры
	coefs = calcDerivativeCoef<double, N>(points); //структура - результат работы функции
	std::cout << "CentralCoef = " << coefs.centralCoef << std::endl; //коэффициент при f(x0)

	std::cout << "OtherCoefs = {";
	for (int(i) = 0; i < N - 1; i++)
	{
		std::cout << coefs.otherCoefs[i] << " "; //массив коэффициентов для узлов
	}
	std::cout << coefs.otherCoefs[N - 1] << "}" << std::endl;


	/////////////////////////точки для графиков
	const int M = 16;
	double x0 = 1;
	std::array <double, M> h;
	for (int(g) = 0; g < M; g++)
	{
		h[g] = pow(10, -g); //массив из шагов
	}

	std::cout << "Шаг и ошибка:" << std::endl;
	for (int(i) = 0; i < M; i++)
	{
		double sum1 = coefs.centralCoef * exp(x0); //A_central*f(x0)
		for (int(g) = 0; g < N; g++)
		{
			sum1 = sum1 + coefs.otherCoefs[g] * exp(x0 + points[g] * h[i]); //D(f)
		}
		std::cout << std::setprecision(M) << h[i] << " " << abs(exp(x0) - sum1 / h[i]) << std::endl; //шаг и ошибка
	}

	return 0;
}
