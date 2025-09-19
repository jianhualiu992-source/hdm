namespace hengduanmian;
public class hdm
{
    //1
    #region   数字转桩号格式 
    /// <summary>
    /// 将数字转换为 K 格式，例如：12345 -> K12+345
    /// </summary>
    /// <param name="n">输入的数字</param>
    /// <returns>格式化后的字符串</returns>
    internal static string Num2K(double meters)
    {
        int km = (int)Math.Floor(meters / 1000);
        double m = meters - km * 1000;
        m = Math.Round(m, 2);

        if (m == (int)m)
        {
            // 整数米，如 2.00 → k1+002
            return $"K{km}+{(int)m:D3}";
        }
        else
        {
            // 小数米，如 2.35 → k0+002.35
            return $"K{km}+{m:000.00}";
        }
    }
    #endregion



    //2
    #region 格式为 DD.MMSS转弧度
    /// <summary>
    /// 将度分秒（DMS）格式转换为弧度
    /// 假设输入的 dms 是一个浮点数，格式为 DD.MMSS
    /// </summary>
    /// <param name="dms">度分秒格式的数值</param>
    /// <returns>对应的弧度值</returns>
    internal static double DmsToRadians(double dms)
    {
        int degrees = (int)dms;
        double fractional = dms - degrees;
        int minutes = (int)(fractional * 100);
        double seconds = (fractional * 100 - minutes) * 100;

        double totalDegrees = degrees + minutes / 60.0 + seconds / 3600.0;
        return totalDegrees * Math.PI / 180;
    }
    #endregion

    //3
    #region  距离和方位角
    /// <summary>
    /// 计算两点之间的极坐标距离和角度（弧度）
    /// </summary>
    /// <param name="x0">起点 X 坐标</param>
    /// <param name="y0">起点 Y 坐标</param>
    /// <param name="x1">终点 X 坐标</param>
    /// <param name="y1">终点 Y 坐标</param>
    /// <returns>包含两个元素的数组：[距离(cd), 角度(hd]) ]</returns>
    internal static double[] Fwj(double x0, double y0, double x1, double y1)
    {
        double x = x1 - x0;
        double y = y1 - y0;

        double cd = Math.Sqrt(x * x + y * y); // 距离
        double hd = Math.Atan2(y, x);         // 角度（弧度）

        if (hd < 0)
        {
            hd += 2 * Math.PI; // 确保角度在 [0, 2π) 范围内
        }

        return [cd, hd];
    }
    #endregion

    //4
    #region 由桩号计算xy
    /// <summary>
    /// k2xy
    /// </summary>
    /// <param name="xyk">起始弧长</param>
    /// <param name="xyx">起始X坐标</param>
    /// <param name="xyy">起始Y坐标</param>
    /// <param name="xyhd">起始方位角（弧度）</param>
    /// <param name="xycd">当前圆弧半径</param>
    /// <param name="xyqdr">起始圆弧半径</param>
    /// <param name="xyzdr">结束圆弧半径</param>
    /// <param name="xyzy">方向系数</param>
    /// <param name="jsk">目标弧长</param>
    /// <param name="jsb">偏移距离</param>
    /// <param name="jd">偏移角度（度）</param>
    /// <returns>包含目标X, 目标Y, 目标方位角（弧度）的对象，同时附加属性 x, y, fwj（度）, rad（弧度）</returns>
    internal static double[] Zs(double xyk, double xyx, double xyy, double xyhd, double xycd, double xyqdr, double xyzdr, double xyzy, double jsk, double jsb, double jd)
    {
        // 第一个条件分支  圆弧
        if (Math.Abs(xyqdr - xyzdr) < 0.01 && xyqdr != 0)
        {
            double centerX = xyx + xyqdr * Math.Cos(xyhd + xyzy * Math.PI / 2);
            double centerY = xyy + xyqdr * Math.Sin(xyhd + xyzy * Math.PI / 2);
            double deltaArcLength = jsk - xyk;
            double deltaAngle = deltaArcLength / xyqdr * xyzy;
            double targetAzimuth = xyhd + deltaAngle;
            if (targetAzimuth < 0)
                targetAzimuth += 2 * Math.PI;
            double angle = xyhd + xyzy * Math.PI / 2 + deltaAngle;
            double targetX = centerX - xyqdr * Math.Cos(angle) + jsb * Math.Cos(angle - xyzy * Math.PI / 2 + jd * Math.PI / 180);
            double targetY = centerY - xyqdr * Math.Sin(angle) + jsb * Math.Sin(angle - xyzy * Math.PI / 2 + jd * Math.PI / 180);

            return [targetX, targetY, targetAzimuth];
        }

        // 第二个条件分支 直线
        if (xyqdr < 0.01 && xyzdr < 0.01 && xyzy < 0.01)
        {
            double targetX = xyx + (jsk - xyk) * Math.Cos(xyhd) + jsb * Math.Cos(xyhd + jd * Math.PI / 180);
            double targetY = xyy + (jsk - xyk) * Math.Sin(xyhd) + jsb * Math.Sin(xyhd + jd * Math.PI / 180);

            return [targetX, targetY, xyhd];
        }

        // 第三个条件分支（默认情况） 缓和曲线
        if (xyqdr < 0.001) xyqdr = 99999999;
        if (xyzdr < 0.001) xyzdr = 99999999;

        double f0 = xyhd;
        double q = xyzy;
        double c = 1 / xyqdr;
        double d = (xyqdr - xyzdr) / (2 * xycd * xyqdr * xyzdr);

        // 初始化 rr 和 vv 数组，索引从1开始，因此大小为5（索引1-4）
        double[] rr = new double[5];
        double[] vv = new double[5];

        rr[1] = 0.1739274226;
        rr[2] = 0.3260725774;
        rr[3] = rr[2];
        rr[4] = rr[1];

        vv[1] = 0.0694318442;
        vv[2] = 0.3300094782;
        vv[3] = 1 - vv[2];
        vv[4] = 1 - vv[1];

        double w = jsk - xyk;
        double xs = 0;
        double ys = 0;

        for (int i = 1; i < 5; i++)
        {
            double ff = f0 + q * vv[i] * w * (c + vv[i] * w * d);
            xs += rr[i] * Math.Cos(ff);
            ys += rr[i] * Math.Sin(ff);
        }

        double fhz3 = f0 + q * w * (c + w * d);
        if (fhz3 < 0)
            fhz3 += 2 * Math.PI;
        if (fhz3 >= 2 * Math.PI)
            fhz3 -= 2 * Math.PI;

        double fhz1 = xyx + w * xs + jsb * Math.Cos(fhz3 + jd * Math.PI / 180);
        double fhz2 = xyy + w * ys + jsb * Math.Sin(fhz3 + jd * Math.PI / 180);
        return [fhz1, fhz2, fhz3];
    }
    #endregion

    //5
    #region 单条路线单个坐标计算
    /// <summary>
    /// 计算指定里程处的坐标和方位角
    /// </summary>
    /// <param name="pqx">路径点数组，每个元素是一个包含多个属性的数组</param>
    /// <param name="k">目标里程</param>
    /// <param name="b">宽度参数</param>
    /// <param name="z">角度参数</param>
    /// <returns>包含目标X坐标、Y坐标和方位角的数组</returns>
    public static double[] Dantiaoxianludange(double[][] pqx, double k, double b, double z)
    {
        int hang = pqx.Length - 1;
        for (int i = 0; i < pqx.Length; i++)
        {
            double dtk = pqx[i][0];
            double dtx = pqx[i][1];
            double dty = pqx[i][2];
            double dtfwj = pqx[i][3];
            double dtcd = pqx[i][4];
            double dtr1 = pqx[i][5];
            double dtr2 = pqx[i][6];
            double dtzy = pqx[i][7];
            if (k >= dtk && k <= dtk + dtcd)
            {
                double hudu = DmsToRadians(dtfwj);
                double[] jsxy1 = Zs(dtk, dtx, dty, hudu, dtcd, dtr1, dtr2, dtzy, k, b, z);
                return [Math.Round(jsxy1[0], 3), Math.Round(jsxy1[1], 3), jsxy1[2]];
            }
        }

        // 如果没有找到对应的里程段，返回起点的坐标和方位角
        return [0, 0, 0];
    }
    #endregion

    //6
    #region xy=>k
    public static double[] Fs(double[][] pqx, double fsx, double fsy)
    {
        dynamic jljd = Fwj(pqx[0][1], pqx[0][2], fsx, fsy);
        double k = pqx[0][0];
        double hudu = DmsToRadians(pqx[0][3]);
        double cz = jljd[0] * Math.Cos(jljd[1] - hudu);
        double pj = jljd[0] * Math.Sin(jljd[1] - hudu);
        int hang = pqx.Length - 1;
        double qdlc = pqx[0][0];
        double zdlc = pqx[hang][0] + pqx[hang][4];
        int jisuancishu = 0;

        while (Math.Abs(cz) > 0.01)
        {
            k = k + cz;
            jisuancishu += 1;

            if (k < qdlc)
            {
                return [-1, -1];
            }

            if (k > zdlc)
            {
                return [-1, -1];
            }

            if (jisuancishu > 15)
            {
                return [-1, -1];
            }

            dynamic xy = Dantiaoxianludange(pqx, k, 0, 0);
            jljd = Fwj((double)xy[0], (double)xy[1], fsx, fsy);
            cz = jljd[0] * Math.Cos(jljd[1] - xy[2]);
            pj = jljd[0] * Math.Sin(jljd[1] - xy[2]);
        }

        //fsjieguo.k = k;
        //fsjieguo.b = pj;
        //fsjieguo.cs = jisuancishu;
        return [Math.Round(k, 3), Math.Round(pj, 3)];
    }
    #endregion

    #region gao
    internal static double Gaocheng(double bpdlc, double bpdgc, double r, double qp, double hp, double t, double k)
    {
        double f = qp - hp;
        r = r * Math.Abs(f) / f;
        double x;
        if (k <= bpdlc - t)
        {
            x = 0;
        }
        else if (k >= bpdlc + t)
        {
            x = 0;
            qp = hp;
        }
        else
        {
            x = k - bpdlc + t;
        }

        return (bpdgc - (bpdlc - k) * qp - Math.Pow(x, 2) / 2 / r);
    }
    #endregion

    //7
    #region h
    public static double h(double k, double[][] sqxb)
    {
        double hp = 0;
        for (int i = 1; i < sqxb.Length - 1; i++)
        {
            double r = sqxb[i][2];
            if (r < 0.001)
                r = 0.001;
            double qp = (sqxb[i][1] - sqxb[i - 1][1]) / (sqxb[i][0] - sqxb[i - 1][0]);
            hp = (sqxb[i + 1][1] - sqxb[i][1]) / (sqxb[i + 1][0] - sqxb[i][0]);
            double f = qp - hp;
            double t = r * Math.Abs(f) / 2;
            if (k <= sqxb[i][0] + t)
                return Math.Round(Gaocheng(sqxb[i][0], sqxb[i][1], r, qp, hp, t, k), 3);
        }
        //the last
        if (k <= sqxb[sqxb.Length - 1][0])
        {
            return Math.Round(sqxb[sqxb.Length - 1][1] + (k - sqxb[sqxb.Length - 1][0]) * hp, 3);
        }
        return -1;
    }
    #endregion


    //polyline offset
    #region

    public static List<(double X, double Y)> OffsetPolyline(double[] xCoords, double[] yCoords, double d)
    {
        int n = xCoords.Length;
        if (n < 2 || yCoords.Length != n)
            throw new ArgumentException("输入坐标数组长度必须相等且至少包含2个点");

        var result = new List<(double X, double Y)>();

        for (int i = 0; i < n - 1; i++)
        {
            (double X, double Y) p1 = (xCoords[i], yCoords[i]);
            (double X, double Y) p2 = (xCoords[i + 1], yCoords[i + 1]);

            double dx = p2.X - p1.X;
            double dy = p2.Y - p1.Y;
            double len = Math.Sqrt(dx * dx + dy * dy);

            if (len < 1e-10)
            {
                if (i == 0) result.Add(p1);
                if (i == n - 2) result.Add(p2);
                continue;
            }

            double ux = dx / len;
            double uy = dy / len;
            double nx = -uy;
            double ny = ux;

            (double X, double Y) offsetStart = (p1.X + nx * d, p1.Y + ny * d);
            (double X, double Y) offsetEnd = (p2.X + nx * d, p2.Y + ny * d);

            if (i == 0)
            {
                result.Add(offsetStart);
            }

            if (i < n - 2)
            {
                (double X, double Y) p3 = (xCoords[i + 2], yCoords[i + 2]);
                double dx2 = p3.X - p2.X;
                double dy2 = p3.Y - p2.Y;
                double len2 = Math.Sqrt(dx2 * dx2 + dy2 * dy2);

                if (len2 < 1e-10)
                {
                    result.Add(offsetEnd);
                    continue;
                }

                double ux2 = dx2 / len2;
                double uy2 = dy2 / len2;
                double nx2 = -uy2;
                double ny2 = ux2;

                (double X, double Y) offset2Start = (p2.X + nx2 * d, p2.Y + ny2 * d);
                (double X, double Y) offset2End = (p3.X + nx2 * d, p3.Y + ny2 * d);

                double x1 = offsetStart.X, y1 = offsetStart.Y;
                double x2 = offsetEnd.X, y2 = offsetEnd.Y;
                double x3 = offset2Start.X, y3 = offset2Start.Y;
                double x4 = offset2End.X, y4 = offset2End.Y;

                double denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);

                if (Math.Abs(denom) > 1e-10)
                {
                    double t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denom;
                    double x = x1 + t * (x2 - x1);
                    double y = y1 + t * (y2 - y1);
                    result.Add((x, y));
                }
                else
                {
                    result.Add(offsetEnd);
                }
            }
            else
            {
                result.Add(offsetEnd);
            }
        }

        return result;
    }

    #endregion


    //fromXgetY
    #region
    /// <summary>
    /// 根据给定的x坐标数组和y坐标数组，计算指定x值对应的y值。
    /// 当x值相等时，返回较大的y值。
    /// </summary>
    /// <param name="x">已知点的x坐标数组（必须升序排列）</param>
    /// <param name="y">已知点的y坐标数组</param>
    /// <param name="targetX">要计算的目标x值</param>
    /// <param name="tolerance">判断x值相等的容差（默认1e-6）</param>
    /// <returns>目标x值对应的y值</returns>
    internal static double fromXgetY(double[] x, double[] y, double targetX)
    {
        double tolerance = 1e-6;
        //  使用二分查找确定targetX所在区间[6](@ref)
        int index = Array.BinarySearch(x, targetX);
        // 如果恰好找到目标x值，直接返回对应的y值
        if (index >= 0)
        {
            return y[index];
        }
        int rightIndex = ~index;
        int leftIndex = rightIndex - 1;
        double xLeft = x[leftIndex];
        double xRight = x[rightIndex];
        if (Math.Abs(xRight - xLeft) < tolerance)
        {
            return Math.Max(y[leftIndex], y[rightIndex]);
        }
        double yLeft = y[leftIndex];
        double yRight = y[rightIndex];
        return yLeft + (yRight - yLeft) * (targetX - xLeft) / (xRight - xLeft);
    }
    #endregion


    #region
    internal static double[] cutAndFillArea(double[] xsA, double[] ysA, double[] xsB, double[] ysB, double epsilon, double jgx = 0)
    {
        //dynamic acad = AcadHelper.GetActiveAutoCAD();
        //acad.Visible = true;
        //dynamic layerFILL = acad.ActiveDocument.Layers.Add("layerFILL");
        //dynamic layerCUT = acad.ActiveDocument.Layers.Add("layerCUT");

        #region
        //extend sj
        double[] left = ExtendLine(xsA[0], ysA[0], xsA[1], ysA[1], -epsilon);
        xsA[0] = left[0]; ysA[0] = left[1];
        //extend dm
        left = ExtendLine(xsB[0], ysB[0], xsB[1], ysB[1], -epsilon);
        xsB[0] = left[0]; ysB[0] = left[1];
        // right extend sj
        //left = ExtendLine(xsA[1], ysA[^1], xsA[^2], ysA[^2], -epsilon);
        //xsA[^1] = left[0]; ysA[^1] = left[1];
        ////right extend dm
        //left = ExtendLine(xsB[^1], ysB[^1], xsB[^2], ysB[^2], -epsilon);
        //xsB[^1] = left[0]; ysB[^1] = left[1];
        left = ExtendLine(xsA[xsA.Length-1], ysA[xsA.Length-1], xsA[xsA.Length-2], ysA[xsA.Length - 2], -epsilon);
        xsA[xsA.Length-1] = left[0]; ysA[xsA.Length - 1] = left[1];
        //right extend dm
        left = ExtendLine(xsB[xsA.Length - 1], ysB[xsA.Length - 1], xsB[xsA.Length - 2], ysB[xsA.Length - 2], -epsilon);
        xsB[xsA.Length - 1] = left[0]; ysB[xsA.Length - 1] = left[1];
        #endregion

        //jgx = Math.Abs(xsB[^1] - xsB[0]);
        //jd list
        List<(double x, double y, int sji, int dmj)> xys = [];
        //area 
        double fill = 0;
        double cut = 0;
        for (int i = 0; i < xsA.Length - 1; i++)
        {
            double x1 = xsA[i];
            double y1 = ysA[i];
            double x2 = xsA[i + 1];
            double y2 = ysA[i + 1];
            for (int j = 0; j < xsB.Length - 1; j++)
            {
                double x3 = xsB[j];
                double y3 = ysB[j];
                double x4 = xsB[j + 1];
                double y4 = ysB[j + 1];
                double[] jd;
                jd = Intersections_Segment(x1, y1, x2, y2, x3, y3, x4, y4);
                if (jd.Length > 0)
                {
                    xys.Add((jd[0], jd[1], i, j));
                }
            }
        }
        //loop xys
        for (int i = 0; i < xys.Count - 1; i++)
        {
            var xy0 = xys[i];
            var xy1 = xys[i + 1];
            int i0 = xy0.sji;
            int i1 = xy1.sji;
            int j0 = xy0.dmj;
            int j1 = xy1.dmj;
            var pts = new List<(double x, double y)>();
            pts.Add((xy0.x, xy0.y));
            //add sj
            for (int j = i0 + 1; j <= i1; j++)
            {
                pts.Add((xsA[j], ysA[j]));
            }
            pts.Add((xy1.x, xy1.y));
            //add dm
            for (int j = j1; j > j0; j--)
            {
                pts.Add((xsB[j], ysB[j]));
            }

            //close if no also ok
            pts.Add((xy0.x, xy0.y));
            //area 
            double area = 0;
            for (int k = 0; k < pts.Count - 1; k++)
                area += pts[k].x * pts[k + 1].y - pts[k].y * pts[k + 1].x;
            double noneed = area < 0 ? fill += Math.Abs(area) / 2.0 : cut += area / 2.0;

            //cad add AcadLWPolyline
            //try
            //{
            //    List<double> lsdmx = new List<double>();
            //    foreach (var pt in pts)
            //    {
            //        lsdmx.Add(pt.x+jgx);
            //        Console.WriteLine(jgx);
            //        lsdmx.Add(pt.y);
            //    }
            //    var dmx = lsdmx.ToArray();
            //    if (area < 0)
            //        acad.ActiveDocument.ActiveLayer = layerFILL;
            //    else
            //        acad.ActiveDocument.ActiveLayer = layerCUT;

            //    dynamic cadline = acad.ActiveDocument.ModelSpace.AddLightWeightPolyline(dmx);
            //    cadline.Closed = 1;
            //}
            //catch
            //{
            //    Console.WriteLine("autocad未运行...");
            //}
        }
        return [Math.Round(fill, 4), Math.Round(cut, 4)];
    }
    #endregion


    //线段
    #region
    /// <summary>
    /// 两线段交点,如果没有则为[]空
    /// </summary>
    /// <param name="x1"></param>
    /// <param name="y1"></param>
    /// <param name="x2"></param>
    /// <param name="y2"></param>
    /// <param name="x3"></param>
    /// <param name="y3"></param>
    /// <param name="x4"></param>
    /// <param name="y4"></param>
    /// <returns></returns>
    [ExcelFunction(Description = "线段延长", Category = "横断面相关")]
    public static double[] Intersections_Segment(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4)
    {
        double d1 = (x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3); // 点(x1,y1)相对于线段(x3,y3)-(x4,y4)
        double d2 = (x4 - x3) * (y2 - y3) - (y4 - y3) * (x2 - x3); // 点(x2,y2)相对于线段(x3,y3)-(x4,y4)
        double d3 = (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1); // 点(x3,y3)相对于线段(x1,y1)-(x2,y2)
        double d4 = (x2 - x1) * (y4 - y1) - (y2 - y1) * (x4 - x1); // 点(x4,y4)相对于线段(x1,y1)-(x2,y2)

        bool isStrictlyIntersecting = (d1 * d2 <= 0) && (d3 * d4 <= 0);

        if (!isStrictlyIntersecting)
            return [];


        double a1 = y2 - y1;
        double b1 = x1 - x2;
        double c1 = x2 * y1 - x1 * y2;
        double a2 = y4 - y3;
        double b2 = x3 - x4;
        double c2 = x4 * y3 - x3 * y4;
        double denominator = a1 * b2 - a2 * b1;
        double px = (b1 * c2 - b2 * c1) / denominator;
        double py = (a2 * c1 - a1 * c2) / denominator;
        return [px, py];
    }
    #endregion


    //延长
    #region
    [ExcelFunction(Description = "线段延长", Category = "横断面相关")]
    public static double[] ExtendLine(double x1, double y1, double x2, double y2, double distance)
    {
        double dx = x2 - x1;
        double dy = y2 - y1;
        double length = Math.Sqrt(dx * dx + dy * dy);
        if (length == 0)
        {
            return [x1, y1];
        }
        double ux = dx / length;
        double uy = dy / length;
        double newX; double newY;
        if (distance < 0)
        {
            newX = x1 + distance * ux;
            newY = y1 + distance * uy;
        }
        else
        {
            newX = x2 + distance * ux;
            newY = y2 + distance * uy;
        }
        return [newX, newY];
    }
    #endregion

    //slope
    #region
    public static double[] slope(double k, double[][] slope)
    {
        for (var i = 0; i < slope.Length - 1; i++)
        {
            if (k >= slope[i][0] && k <= slope[i + 1][0])
            {
                double L = slope[i][1]
                        + (slope[i + 1][1] - slope[i][1]) / (slope[i + 1][0] - slope[i][0]) * (k - slope[i][0]);
                double R = slope[i][2]
                        + (slope[i + 1][2] - slope[i][2]) / (slope[i + 1][0] - slope[i][0]) * (k - slope[i][0]);
                return [Math.Round(L, 3), Math.Round(R, 3)];
            }
        }
        return [0, 0];
    }
    #endregion

    //根据给定的x坐标数组和y坐标数组，计算指定x值对应的y值。
    #region
    /// <summary>
    /// 根据给定的x坐标数组和y坐标数组，计算指定x值对应的y值。
    /// 当x值相等时，返回较大的y值。
    /// </summary>
    /// <param name="x">已知点的x坐标数组（必须升序排列）</param>
    /// <param name="y">已知点的y坐标数组</param>
    /// <param name="targetX">要计算的目标x值</param>
    /// <param name="tolerance">判断x值相等的容差（默认1e-6）</param>
    /// <returns>目标x值对应的y值</returns>
    internal static double Interpolate(double[] x, double[] y, double targetX)
    {
        double tolerance = 1e-6;
        //  使用二分查找确定targetX所在区间[6](@ref)
        int index = Array.BinarySearch(x, targetX);
        // 如果恰好找到目标x值，直接返回对应的y值
        if (index >= 0)
        {
            return y[index];
        }
        int rightIndex = ~index;
        int leftIndex = rightIndex - 1;
        double xLeft = x[leftIndex];
        double xRight = x[rightIndex];
        if (Math.Abs(xRight - xLeft) < tolerance)
        {
            return Math.Max(y[leftIndex], y[rightIndex]);
        }
        double yLeft = y[leftIndex];
        double yRight = y[rightIndex];
        return yLeft + (yRight - yLeft) * (targetX - xLeft) / (xRight - xLeft);
    }
    #endregion


    // 折线交点
    #region
    [ExcelFunction(Description = "折线交点", Category = "横断面相关")]
    public static double[,] Intersections_polyline(double[] xsA, double[] ysA, double[] xsB, double[] ysB, double epsilon)
    {
        #region
        //extend sj
        if (Math.Abs(epsilon) > 0.001)
        {
            double[] left = ExtendLine(xsA[0], ysA[0], xsA[1], ysA[1], -epsilon);
            xsA[0] = left[0]; ysA[0] = left[1];
            //extend dm
            left = ExtendLine(xsB[0], ysB[0], xsB[1], ysB[1], -epsilon);
            xsB[0] = left[0]; ysB[0] = left[1];
            // right extend sj
            left = ExtendLine(xsA[xsA.Length - 1], ysA[xsA.Length - 1], xsA[xsA.Length - 2], ysA[xsA.Length - 2], -epsilon);
            xsA[xsA.Length - 1] = left[0]; ysA[xsA.Length - 1] = left[1];
            //right extend dm
            left = ExtendLine(xsB[xsA.Length - 1], ysB[xsA.Length - 1], xsB[xsA.Length - 2], ysB[xsA.Length - 2], -epsilon);
            xsB[xsA.Length - 1] = left[0]; ysB[xsA.Length - 1] = left[1];
        }
        #endregion
        //jd list
        List<(double x, double y, int sji, int dmj)> xys = [];
        for (int i = 0; i < xsA.Length - 1; i++)
        {
            double x1 = xsA[i];
            double y1 = ysA[i];
            double x2 = xsA[i + 1];
            double y2 = ysA[i + 1];
            for (int j = 0; j < xsB.Length - 1; j++)
            {
                double x3 = xsB[j];
                double y3 = ysB[j];
                double x4 = xsB[j + 1];
                double y4 = ysB[j + 1];
                double[] jd;
                jd = Intersections_Segment(x1, y1, x2, y2, x3, y3, x4, y4);
                if (jd.Length > 0)
                {
                    xys.Add((jd[0], jd[1], i, j));
                }
            }
        }
        double[,] xy = new double[xys.Count, 4];
        for (int i = 0; i < xys.Count; i++)
        {
            xy[i, 0] = xys[i].x;
            xy[i, 1] = xys[i].y;
            xy[i, 2] = xys[i].sji;
            xy[i, 3] = xys[i].dmj;
        }
        return xy;
    }
    #endregion


    //直线交点
    [ExcelFunction(Description = "直线交点", Category = "横断面相关")]
    public static double[] Intersections_line(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4)
    {
        double a1 = y2 - y1;
        double b1 = x1 - x2;
        double c1 = x2 * y1 - x1 * y2;
        double a2 = y4 - y3;
        double b2 = x3 - x4;
        double c2 = x4 * y3 - x3 * y4;
        double denominator = a1 * b2 - a2 * b1;
        if (Math.Abs(denominator) < 1e-10)
            return [0,0]; 
        double px = (b1 * c2 - b2 * c1) / denominator;
        double py = (a2 * c1 - a1 * c2) / denominator;
        return [px, py];
    }


}
