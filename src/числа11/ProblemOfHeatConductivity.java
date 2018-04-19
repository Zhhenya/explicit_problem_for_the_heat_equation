package числа11;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class ProblemOfHeatConductivity {
    private List<Double> yExplict;
    private List<Double> yObvious;
    private List<Double> prevY;
    private List<Double> realY;
    private List<Double> coeffA;
    private List<Double> coeffB;
    private List<Double> coeffC;
    private List<Double> coeffD;
    private List<Double> x;
    private List<Double> t;
    private List<Double> fi;
    private int countOfPart = 200;
    private int n = countOfPart + 1;
    private double T = 0.02;
    private int numberOfLay = 5;
    private int countOfLay = 0;

    private double p = 20;
    private double a = Math.sqrt(p)/4;
    private double l = Math.PI/2;
    private double alfa1 = 0;
    private double beta1 = 1;
    private double alfa2 = 1;
    private double beta2 = 0;
    private double Mu1 = 0;
    private double Mu2 = 0;
    private double sigma = 1;

    private double h;
    private double tao;

    private void initAll(){

        yExplict = new ArrayList<>();
        yObvious = new ArrayList<>();
        prevY = new ArrayList<>();
        realY = new ArrayList<>();
        coeffA = new ArrayList<>();
        coeffB = new ArrayList<>();
        coeffC = new ArrayList<>();
        coeffD = new ArrayList<>();
        fi = new ArrayList<>();
        x = new ArrayList<>();
        t = new ArrayList<>();
    }

    private double u0(double x){
        return Math.sin(x);
    }
    private double f(double t, double x){
        return (p/16) * Math.exp(p * t/16) * Math.sin(x);
    }
    private double u(double t, double x){
        return ((Math.exp(p * t/16) + Math.exp(-p * t/16))/2) * Math.sin(x);
    }
    private double KSI1 = 0;
    private double KSI2(){
        return alfa2*sigma/delta2();
    }

    private double v1(){
        return Mu1 /beta1;
    }
    private double delta1(){
        return sigma * (alfa1 + h * beta1) + alfa1 * Math.pow(h, 2)/(2 * Math.pow(a, 2) * tao);
    }
    private double delta2(){
        return sigma * (alfa2 + h * beta2) + alfa2 * Math.pow(h, 2)/(2 * Math.pow(a, 2) * tao);
    }
    private double Mu1_(double t){
        return Mu1 + alfa1 * h * f(t, 0)/(2 * Math.pow(a, 2));
    }
    private double Mu2_(double t){
        return Mu2 + alfa2 * h * f(t, l)/(2 * Math.pow(a, 2));
    }
    private double v2(double t, double prevY, double y){
        return 1/delta2() * (h * Mu2_(t) + (1 - sigma)*alfa2 * prevY + y * (alfa2 * Math.pow(h, 2)/(2 * Math.pow(a, 2) * tao) - (1 - sigma) * (alfa2 + beta2 * h)));
    }
    private double fi(double t, double x){
        return f(t + sigma * tao, x);
    }
    private double A(){
        return Math.pow(a, 2) * sigma;
    }
    private double C;
    private double B(){
        return 2 * Math.pow(a, 2) * sigma + Math.pow(h, 2)/tao;
    }
    private double D(double x, double t, double y, double prevY, double nextY, double fi){
        return (2 * Math.pow(a, 2) * (1 - sigma) - Math.pow(h, 2)/tao) * y - Math.pow(a, 2) * (1 - sigma) * (prevY + nextY) - fi * Math.pow(h, 2);
    }

    public void initH(){
        h = l/countOfPart;
        tao = T/numberOfLay;
    }
    public void initStartY(){
        for(int i = 0; i < n; i++) {
            prevY.add(u0(x.get(i)));
            yExplict.add(0.0);
        }
    }

    public void initXAndT(){
        for(int i =  0; i <  n; i++)
            x.add(i * h);
        for(int i = 0; i <= numberOfLay; i++)
            t.add(i * tao);


    }


    public void initCoefForNotAnExplictMethod(){
        coeffA.add(0.0);
    //    coeffB.add(1.0);
        coeffB.add(beta1);
        coeffC.add(KSI1);
        coeffD.add(v1());
        for(int i  = 1; i < n - 1; i++){
            coeffA.add(A());
            coeffB.add(B());
            coeffC.add(C);
            coeffD.add(D(x.get(i), t.get(countOfLay), prevY.get(i), prevY.get(i - 1), prevY.get(i + 1), fi(t.get(countOfLay), x.get(i))));
        }
        coeffA.add(-KSI2());
        coeffB.add(-1.0);
        coeffC.add(0.0);
        coeffD.add(v2(t.get(countOfLay), prevY.get(n - 2), prevY.get(n - 1)));
    }

    public void solutionForAll(){
        C = A();
        initAll();
        initH();
        initXAndT();
        T = 1;
        double r = 0;
        do {
            tao = T / numberOfLay;
            r = Math.pow(a, 2) * tao / Math.pow(h, 2);
            sigma = 0;
            if(r > 0.5)
                T = T /10;
        }while(r > 0.5);
        sigma = 1;
        solutionForImplicitMethod();
        sigma = 0;
        solutionForObvious(r);
        print();
        norm(yExplict, realY, "неявного метода");
        norm(yObvious, realY, "явного метода");

    }

    public void solutionForObvious(double r){

        initStartY();
        countOfLay = 1;
        while(countOfLay < numberOfLay){
            yObvious.clear();
            yObvious.add(v1());
            for(int i = 1; i < n - 1; i++)
                yObvious.add((1 - 2 * r) * prevY.get(i) + r * (prevY.get(i - 1) + prevY.get(i + 1)) + tao * fi(t.get(countOfLay), x.get(i)));
            yObvious.add(v2(KSI2() * yObvious.get(n - 2) + t.get(countOfLay + 1), prevY.get(n - 2), prevY.get(n - 1)));
            countOfLay++;
            Collections.copy(prevY, yObvious);
        }
        getRealY();
        Collections.copy(yObvious, prevY);
    }

    public void solutionForImplicitMethod(){
        countOfLay  =1;
        initStartY();

        while(countOfLay < numberOfLay) {
            initCoefForNotAnExplictMethod();
            Run.initCoefForEquation(coeffA, coeffB, coeffC, coeffD);
            Run.initStartCoeffForMethod(n);
            Run.tridiagonalMatrixSystem(n);
            Collections.copy(prevY, Run.getY());
            Run.clear();
            countOfLay++;
        }
        Collections.copy(yExplict, prevY);
        prevY.clear();

    }

    public void norm(List<Double> y, List<Double> realY, String s){
        double max = 0;
        for(int i = 0; i < n; i++)
            if(Math.abs((y.get(i) - realY.get(i))/y.get(i)) > max)
                max = Math.abs((y.get(i) - realY.get(i))/y.get(i));
        System.out.println("Норма относительной погрешности для " + s + " c шагом " + h + " = " + max);
    }
    public void getRealY(){
        realY.clear();
        for(int  i = 0; i < n; i++)
            realY.add(u(t.get(countOfLay), x.get(i)));
    }

    public void print(){
        System.out.println("\t\t\t\t X \t\t\t\t  Yнеявный \t\t\t\t Yявный\t\t\t\tReal Y");
        for(int  i = 0; i < x.size(); i++)
            System.out.printf("%20.10f %20.10f %20.10f %20.10f \n",x.get(i) , yExplict.get(i), yObvious.get(i),  realY.get(i));
    }
}