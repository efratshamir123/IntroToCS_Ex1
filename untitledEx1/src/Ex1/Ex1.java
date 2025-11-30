package Ex1;

import java.util.Arrays;

/**
 * Introduction to Computer Science 2026, Ariel University,
 * Ex1: arrays, static functions and JUnit
 * https://docs.google.com/document/d/1GcNQht9rsVVSt153Y8pFPqXJVju56CY4/edit?usp=sharing&ouid=113711744349547563645&rtpof=true&sd=true
 *
 * This class represents a set of static methods on a polynomial functions - represented as an array of doubles.
 * The array {0.1, 0, -3, 0.2} represents the following polynomial function: 0.2x^3-3x^2+0.1
 * This is the main Class you should implement (see "add your code below")
 *
 * @author boaz.benmoshe

 */
public class Ex1 {
	/** Epsilon value for numerical computation, it serves as a "close enough" threshold. */
	public static final double EPS = 0.001; // the epsilon to be used for the root approximation.
	/** The zero polynomial function is represented as an array with a single (0) entry. */
	public static final double[] ZERO = {0};
	/**
	 * Computes the f(x) value of the polynomial function at x.
	 * @param poly - polynomial function
	 * @param x
	 * @return f(x) - the polynomial function value at x.
	 */
	public static double f(double[] poly, double x) {
		double ans = 0;
		for(int i=0;i<poly.length;i++) {
			double c = Math.pow(x, i);
			ans += c*poly[i];
		}
		return ans;
	}
	/** Given a polynomial function (p), a range [x1,x2] and an epsilon eps.
	 * This function computes an x value (x1<=x<=x2) for which |p(x)| < eps, 
	 * assuming p(x1)*p(x2) <= 0.
	 * This function should be implemented recursively.
	 * @param p - the polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
	 * @return an x value (x1<=x<=x2) for which |p(x)| < eps.
	 */
	public static double root_rec(double[] p, double x1, double x2, double eps) {
		double f1 = f(p,x1);
		double x12 = (x1+x2)/2;
		double f12 = f(p,x12);
		if (Math.abs(f12)<eps) {return x12;}
		if(f12*f1<=0) {return root_rec(p, x1, x12, eps);}
		else {return root_rec(p, x12, x2, eps);}
	}
	/**
	 * This function computes a polynomial representation from a set of 2D points on the polynom.
	 * The solution is based on: //	http://stackoverflow.com/questions/717762/how-to-calculate-the-vertex-of-a-parabola-given-three-points
	 * Note: this function only works for a set of points containing up to 3 points, else returns null.
	 * @param xx
	 * @param yy
	 * @return an array of doubles representing the coefficients of the polynom.
	 */
	public static double[] PolynomFromPoints(double[] xx, double[] yy) {
		double [] ans = null;
		int lx = xx.length;
		int ly = yy.length;
		if(xx!=null && yy!=null && lx==ly && lx>1 && lx<4) {
            /** add you code below

             /////////////////// */
            if(lx == 2) {
                    double a = (yy[1] - yy[0]) / (xx[1] - xx[0]);
                    double b = yy[0] - a * xx[0];
                    ans = new double[]{b, a};
                }

                else if(lx == 3) {
                    double x1 = xx[0], x2 = xx[1], x3 = xx[2];
                    double y1 = yy[0], y2 = yy[1], y3 = yy[2];

                    double d = (x1 - x2)*(x1 - x3)*(x2 - x3);

                    double a = (y1*(x2 - x3) + y2*(x3 - x1) + y3*(x1 - x2)) / d;
                    double b = (y1*(x3*x3 - x2*x2) + y2*(x1*x1 - x3*x3) + y3*(x2*x2 - x1*x1)) / d;
                    double c = (y1*(x2*x3*(x2 - x3)) + y2*(x3*x1*(x3 - x1)) + y3*(x1*x2*(x1 - x2))) / d;

                    ans = new double[]{c, b, a};
                }
        }
        return ans;
	}
	/** Two polynomials functions are equal if and only if they have the same values f(x) for n+1 values of x,
	 * where n is the max degree (over p1, p2) - up to an epsilon (aka EPS) value.
	 * @param p1 first polynomial function
	 * @param p2 second polynomial function
	 * @return true iff p1 represents the same polynomial function as p2.
	 */
	public static boolean equals(double[] p1, double[] p2) {
		boolean ans = true;
        /** add you code below
         /////////////////// */
        int n = Math.max(p1.length, p2.length);
        for(int i=0; i<=n; i++){
            double x = i;
            double f1 = f(p1, x);
            double f2 = f(p2, x);
            if(Math.abs(f1 - f2) > EPS) ans = false;
        }

        return ans;
	}

	/** 
	 * Computes a String representing the polynomial function.
	 * For example the array {2,0,3.1,-1.2} will be presented as the following String  "-1.2x^3 +3.1x^2 +2.0"
	 * @param poly the polynomial function represented as an array of doubles
	 * @return String representing the polynomial function:
	 */
	public static String poly(double[] poly) {
		String ans = "";
		if(poly.length==0) {ans="0";}
		else {
            /** add you code below
             /////////////////// */
            for(int i = poly.length - 1; i >= 0; i--) {
                double c = poly[i];
                if(Math.abs(c) > EPS) {
                    if(!ans.equals("") && c > 0) ans += " +";
                    if(c < 0) ans += " ";

                    ans += c;

                    if(i == 1) ans += "x";
                    if(i > 1) ans += "x^" + i;
                }
            }
            if(ans.equals("")) ans = "0";
		}
		return ans;
	}
	/**
	 * Given two polynomial functions (p1,p2), a range [x1,x2] and an epsilon eps. This function computes an x value (x1<=x<=x2)
	 * for which |p1(x) -p2(x)| < eps, assuming (p1(x1)-p2(x1)) * (p1(x2)-p2(x2)) <= 0.
	 * @param p1 - first polynomial function
	 * @param p2 - second polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
	 * @return an x value (x1<=x<=x2) for which |p1(x) - p2(x)| < eps.
	 */
	public static double sameValue(double[] p1, double[] p2, double x1, double x2, double eps) {
		double ans = x1;
        /** add you code below
         /////////////////// */
        double left = x1;
        double right = x2;
        double fLeft = f(p1, left) - f(p2, left);

        for (int i = 0; i < 100; i++) {
            double mid = (left + right) / 2;
            double fMid = f(p1, mid) - f(p2, mid);
            ans = mid;

            if (Math.abs(fMid) < eps) {
                break;
            }

            if (fLeft * fMid <= 0) {
                right = mid;
            } else {
                left = mid;
                fLeft = fMid;
            }
        }

		return ans;
	}
	/**
	 * Given a polynomial function (p), a range [x1,x2] and an integer with the number (n) of sample points.
	 * This function computes an approximation of the length of the function between f(x1) and f(x2) 
	 * using n inner sample points and computing the segment-path between them.
	 * assuming x1 < x2. 
	 * This function should be implemented iteratively (none recursive).
	 * @param p - the polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param numberOfSegments - (A positive integer value (1,2,...).
	 * @return the length approximation of the function between f(x1) and f(x2).
	 */
	public static double length(double[] p, double x1, double x2, int numberOfSegments) {
		double ans = x1;
        /** add you code below
         /////////////////// */
        double dx = (x2 - x1) / numberOfSegments;          // אורך של כל קטע
        double x = x1;                                     // הנקודה הקודמת מתחילה ב x1
        double y = f(p, x);                               //ערך הפונקציה בנקודה הקודמת

        ans = 0;                                          // מתחילים את סכום האורך מ-0

        for (int i = 1; i <= numberOfSegments; i = i +1){ // לולאה שעוברת על כל הקטעים

            double xNew = x1 + i * dx;                    // מחשב את הנקודה הבאה
            double yNew = f(p, xNew);                     // מחשב את ערך הפולינום בנקודה החדשה

            double dxDistance = xNew - x;                 // הפרש ב-x בין נקודה חדשה ונקודה קודמת
            double dyDistance = yNew - y;                 // הפרש ב-y בין נקודה חדשה ונקודה קודמת

            double distance = Math.sqrt(dxDistance*dxDistance + dyDistance*dyDistance); // נוסתחת distance
            ans = ans + distance;                         // מוסיפים את אורך הקטע לסכום הכללי

            x = xNew;                                     // מעדכנים את הנקודה הקודמת
            y = yNew;
        }
		return ans;
	}
	
	/**
	 * Given two polynomial functions (p1,p2), a range [x1,x2] and an integer representing the number of Trapezoids between the functions (number of samples in on each polynom).
	 * This function computes an approximation of the area between the polynomial functions within the x-range.
	 * The area is computed using Riemann's like integral (https://en.wikipedia.org/wiki/Riemann_integral)
	 * @param p1 - first polynomial function
	 * @param p2 - second polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param numberOfTrapezoid - a natural number representing the number of Trapezoids between x1 and x2.
	 * @return the approximated area between the two polynomial functions within the [x1,x2] range.
	 */
	public static double area(double[] p1,double[]p2, double x1, double x2, int numberOfTrapezoid) {
		double ans = 0;
        /** add you code below
         /////////////////// */
        if (numberOfTrapezoid <= 0) return 0;

        double dx = (x2 - x1) / numberOfTrapezoid;
        double totalArea = 0;

        double prevX = x1;
        double prevY1 = f(p1, prevX);
        double prevY2 = f(p2, prevX);

        for (int i = 1; i <= numberOfTrapezoid; i++) {
            double currX = x1 + i * dx;
            double currY1 = f(p1, currX);
            double currY2 = f(p2, currX);

            if ((prevY1 - prevY2) * (currY1 - currY2) < 0) {

                double crossX = sameValue(p1, p2, prevX, currX, EPS);
                double crossY1 = f(p1, crossX);
                double crossY2 = f(p2, crossX);

                totalArea += 0.5 * (Math.abs(prevY1 - prevY2) + Math.abs(crossY1 - crossY2)) * (crossX - prevX);

                prevX = crossX;
                prevY1 = crossY1;
                prevY2 = crossY2;
            }

            totalArea += 0.5 * (Math.abs(prevY1 - prevY2) + Math.abs(currY1 - currY2)) * (currX - prevX);

            prevX = currX;
            prevY1 = currY1;
            prevY2 = currY2;
        }

        ans = totalArea;

        return ans;
	}
	/**
	 * This function computes the array representation of a polynomial function from a String
	 * representation. Note:given a polynomial function represented as a double array,
	 * getPolynomFromString(poly(p)) should return an array equals to p.
	 * 
	 * @param p - a String representing polynomial function.
	 * @return
	 */
	public static double[] getPolynomFromString(String p) {
		double [] ans = ZERO;//  -1.0x^2 +3.0x +2.0
        /** add you code below
         /////////////////// */
        p = p.replace(" ", "");
        p = p.replace("-", "+-");
        String[] terms = p.split("\\+");

        int maxPow = 0;

        for(String t : terms){
            if(t.equals("")) continue;
            if(t.contains("x")){
                if(t.contains("^")){
                    int pow = Integer.parseInt(t.substring(t.indexOf("^") + 1));
                    if(pow > maxPow) maxPow = pow;
                } else {
                    if(maxPow < 1) maxPow = 1;
                }
            }
        }

        double[] poly = new double[maxPow + 1];

        for(String t : terms){
            if(t.equals("")) continue;

            double coef;
            int pow;

            if(t.contains("x")){
                if(t.equals("x")) coef = 1;
                else if(t.equals("-x")) coef = -1;
                else {
                    int idx = t.indexOf("x");
                    coef = Double.parseDouble(t.substring(0, idx));
                }

                if(t.contains("^")){
                    pow = Integer.parseInt(t.substring(t.indexOf("^") + 1));
                } else {
                    pow = 1;
                }
            }
            else {
                coef = Double.parseDouble(t);
                pow = 0;
            }

            poly[pow] = coef;
        }

        ans = poly;

        return ans;
	}
	/**
	 * This function computes the polynomial function which is the sum of two polynomial functions (p1,p2)
	 * @param p1
	 * @param p2
	 * @return
	 */
	public static double[] add(double[] p1, double[] p2) {
		double [] ans = ZERO;//
        /** add you code below
         /////////////////// */
        int len1 = p1.length;                   // אורך הפולינום p1
        int len2 = p2.length;                   // אורך הפולינום p2
        int maxLen = Math.max(len1, len2);      // לוקחים את האורך הגדול יותר מבין שניהם (הפולינום על החזקה הגבוהה)

        ans = new double[maxLen];               //בונים מערך חדש בגודל האורך המקסימלי

        for (int i = 0; i < maxLen; i = i + 1){ //לולאה שעוברת על כל החזקות מ-0 עד maxlen
            double c1 = 0;                      //מגדירים מקדם ל-p1, אם אין ערך המקדם 0
            double c2 = 0;                      // מגדירים מקדם ל-p2, אם אין ערך המקדם 0

            if (i < len1) c1 = p1[i];           //אם קיים ערך ב-p1 ניקח אותו
            if (i < len2) c2 = p2[i];           //אם קיים ערך ב-p2 ניקח אותו

            ans[i] = c1 + c2;                   // מחברים את המקדמים ושמים במערך באינדקס המתאים
        }
        int newLen = maxLen - 1;
        while (newLen > 0 && ans[newLen] == 0) newLen--;

        ans = Arrays.copyOf(ans, newLen + 1);

		return ans;
	}
	/**
	 * This function computes the polynomial function which is the multiplication of two polynoms (p1,p2)
	 * @param p1
	 * @param p2
	 * @return
	 */
	public static double[] mul(double[] p1, double[] p2) {
		double [] ans = ZERO;
        /** add you code below
         /////////////////// */
        int len1 = p1.length;                              //אורך הפולינום p1
        int len2 = p2.length;                              // אורך הפולינום p2

        int newLen = len1 + len2 - 1;                      // אורך הפולינום החדש אחרי הכפל
        ans = new double[newLen];

        for (int i = 0; i < len1; i = i + 1){              // לולאה שרצה על כל האיברים ב-p1
            for (int j = 0; j < len2; j = j + 1){          // לולאה שרצה על כל האיברים ב-p2
                ans[i + j] += (p1[i] * p2[j]);             //מכפילים את המקדמים לפי i,j ושמים לפי החזקה i+j
            }
        }
        int k = newLen - 1;
        while (k > 0 && ans[k] == 0)
            k--;
        ans = Arrays.copyOf(ans, k + 1);
		return ans;
	}
	/**
	 * This function computes the derivative of the p0 polynomial function.
	 * @param po
	 * @return
	 */
	public static double[] derivative (double[] po) {
		double [] ans = ZERO;//
        /** add you code below
         /////////////////// */
        int len = po.length;                //שומרים את אורך הפולינום המקורי
        if (len == 0) return po;            //אם אין פולינום מחזירים כמו שהוא
        if (len == 1) return ZERO;          // אם יש רק מספר אז הנגזרת 0

        ans = new double[len - 1];          // הנגזרת קצרה יותר באיבר אחד

        for (int i = 1; i < len; i = i+1){  // עוברים על כל החזקות מ-1 כי לחזקה 0 אין נגזרת
            double coef = po[i];            // לוקחים את המקדם של החזקה i
            int power = i;                  // החזקה של האיבר זה האינדקס שלו
            ans[i - 1] = coef * power;      // מציבים את המקדם החדש בחזקה
        }
		return ans;
	}
}
