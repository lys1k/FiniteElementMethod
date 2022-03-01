import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.BaseAbstractUnivariateIntegrator;
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator;
import org.apache.commons.math3.linear.*;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Scanner;

public class Mes {
    public static int N;
    public static int x0 = 0;
    public static int xn = 2;
    public static int betta = 1;
    public static int gamma = 0;
    public static int a = -1;
    public static int b = 0;
    public static int c = -1;
    public static RealVector solution;


    static IterativeLegendreGaussIntegrator integrator = new IterativeLegendreGaussIntegrator(10,
            BaseAbstractUnivariateIntegrator.DEFAULT_RELATIVE_ACCURACY,
            BaseAbstractUnivariateIntegrator.DEFAULT_ABSOLUTE_ACCURACY);


    static double function(double x){
        return Math.sin(x);
    }


    static double base(int i,double x){
        double xi_1 = (double) xn * (i - 1)/N;
        double xi = (double) xn * i/N;
        double xi1 = (double) xn*  (i+1)/N;

        if (x >= xi_1 && x < xi){
            return (x- xi_1) / (xi - xi_1);
        }
        else if (x >= xi && x < xi1){
            return (xi1 - x)/(xi1 - xi);
        }
        return 0.0;
    }


    static double base_der(int i,double x){
        double xi_1 = (double)xn*(i-1)/N;
        double xi = (double)xn*(i)/N;
        double xi1 = (double)xn*(i+1)/N;

        if (x >= xi_1 && x < xi){
            return N / (double) xn;
        }
        else if (x >= xi && x < xi1){
            return -N / (double) xn;
        }
        return 0.0;
    }


    static double calculate_B(int i, int j){
        UnivariateFunction f1 = x -> - a * base_der(i, x) * base_der(j, x);
        UnivariateFunction f2 = x -> b *  base_der(i, x) * base(j,x);
        UnivariateFunction f3 = x -> c * base(i, x) * base(j, x);

        return integrator.integrate(Integer.MAX_VALUE, f1, x0, xn)
                + integrator.integrate(Integer.MAX_VALUE, f2, x0, xn)
                + integrator.integrate(Integer.MAX_VALUE, f3, x0, xn)
                - betta * (base(i,xn)) * base(j, xn);
    }


    static double calculate_L(int i){
        UnivariateFunction f = x -> function(x) * base(i, x);

        return integrator.integrate(Integer.MAX_VALUE, f, x0, xn) - gamma * base(i, xn);
    }


    static double u(double x, double[] solution){
        double multi = 0;
        for(int i = 0; i < N + 1 ; i++){
            multi += solution[i] * base(i, x);
        }
        return multi;
    }


    static void find_solution(double[] l_matrix, double[][] b_matrix){
        RealMatrix B = new Array2DRowRealMatrix(b_matrix);
        DecompositionSolver solver = new LUDecomposition(B).getSolver();
        RealVector L = new ArrayRealVector(l_matrix);
        solution = solver.solve(L);
    }


    static void save_file() throws FileNotFoundException {
        double[] sol = solution.toArray();
        double[] x = new double[200];
        double[] y = new double[200];

        for (int i = 0; i < 200; i++) {
            x[i] = (double) xn * i / 200;
            y[i] = u(x[i], sol);
        }

        PrintWriter writer = new PrintWriter("result.txt");
        for(int i = 0; i < 200; i++) {
            writer.println(x[i] + " " + y[i]);
        }
        writer.close();
    }


    public static void main(String[] args) throws FileNotFoundException {
        Scanner scanner = new Scanner(System.in);
        System.out.println("Liczba przedziałów N:");
        N = scanner.nextInt();

        double[][] b_matrix = new double[N+1][N+1];
        double[] l_matrix = new double[N+1];

        for(int i  = 0; i < N +1 ; i++){
            for(int j = 0; j < N +1 ; j++){
                b_matrix[j][i] = calculate_B(j, i);
            }
        }

        for(int i = 0; i < N+1 ; i++){
            l_matrix[i] = calculate_L(i);
        }

        find_solution(l_matrix, b_matrix);
        save_file(); //zapisuje do pliku, aby móc eksportować dane do excela

        System.out.println("Macierz L:");
        System.out.println(Arrays.toString(l_matrix));

        System.out.println("\nMacierz B:");
        for(int i=0; i<= N; i++){
            for(int j=0; j <= N; j++){
                System.out.print(b_matrix[i][j] + " ");
            }
            System.out.println();
        }

        System.out.println("\nRozwiązanie:");
        System.out.println(solution);
    }
}
