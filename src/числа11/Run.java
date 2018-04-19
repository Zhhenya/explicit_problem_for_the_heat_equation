package числа11;



import java.util.ArrayList;
import java.util.List;

public class Run {
    private static List<Double> coeffAlfa = new ArrayList<>();
    private static List<Double> coeffBeta = new ArrayList<>();
    private static List<Double> y = new ArrayList<>();
    private static List<Double> coeffA = new ArrayList<>();
    private static List<Double> coeffB = new ArrayList<>();
    private static List<Double> coeffC = new ArrayList<>();
    private static List<Double> coeffD = new ArrayList<>();


    public static void initCoefForEquation(List<Double> A, List<Double> B, List<Double> C, List<Double> D){
        coeffA = A;
        coeffB = B;
        coeffC = C;
        coeffD = D;

    }

    public static void initStartCoeffForMethod(int m){
        coeffAlfa.add(0.0);
        coeffBeta.add(0.0);
        coeffAlfa.add(coeffC.get(0)/coeffB.get(0));
        coeffBeta.add( -coeffD.get(0)/coeffB.get(0));

        y.add(0.0);
        y.add(0.0);

        for(int i = 2; i <= m; i++){
            coeffAlfa.add(0.0);
            coeffBeta.add(0.0);
            if(i < m)
                y.add(0.0);
        }
    }

    public static void tridiagonalMatrixSystem(int m) {
        //прямой ход
        for (int i = 2; i <= m; i++) {
            coeffAlfa.set(i, coeffC.get(i - 1) / (coeffB.get(i - 1) - coeffAlfa.get(i - 1) * coeffA.get(i - 1)));
            coeffBeta.set(i, (coeffA.get(i - 1) * coeffBeta.get(i - 1) - coeffD.get(i - 1)) / (coeffB.get(i - 1) - coeffA.get(i - 1) * coeffAlfa.get(i - 1)));
        }
        //обратный ход
        y.set(m - 1, coeffBeta.get(m));
        for (int i = m - 2; i >= 0; i--)
            y.set(i, coeffC.get(i) * y.get(i + 1) / (coeffB.get(i) - coeffA.get(i) * coeffAlfa.get(i)) + (coeffA.get(i) * coeffBeta.get(i) - coeffD.get(i)) / (coeffB.get(i) - coeffA.get(i) * coeffAlfa.get(i)));

    }
    public static List<Double> getY() {
        return y;
    }
    public static void clear(){
        coeffAlfa.clear();
        coeffBeta.clear();
        coeffA.clear();
        coeffB.clear();
        coeffC.clear();
        coeffD.clear();
        y.clear();
    }
}
