using System.Runtime.InteropServices;

namespace hengduanmian
{
    public class AcadHelper
    {
        // 导入 COM 相关 API
        [DllImport("ole32.dll")]
        private static extern int CLSIDFromProgID(
            [MarshalAs(UnmanagedType.LPWStr)] string lpszProgID,
            out Guid pclsid);


        [DllImport("oleaut32.dll")]
        private static extern int GetActiveObject(
            [MarshalAs(UnmanagedType.LPStruct)] Guid rclsid,
            IntPtr pvReserved,
            [MarshalAs(UnmanagedType.Interface)] out object ppunk);


        // 获取正在运行的 AutoCAD 实例
        public static dynamic GetActiveAutoCAD(string progId = "AutoCAD.Application")
        {
            Guid clsid;
            CLSIDFromProgID(progId, out clsid);  // 获取 AutoCAD 的 CLSID
            GetActiveObject(clsid, IntPtr.Zero, out var obj);  // 获取活动对象
            return obj;  // 返回 AutoCAD.Application 对象
        }
    }

}
