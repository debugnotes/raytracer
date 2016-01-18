using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using MatrixUtils;

namespace Raytracer
{
    public partial class MainForm : Form
    {
        public MainForm()
        {
            InitializeComponent();
        }

        private void MainForm_Load(object sender, EventArgs e)
        {
            //MatrixUtils.UnitTesting ut = new MatrixUtils.UnitTesting();
            //ut.AllTests();

            Raytracer rt = new Raytracer();

            //rt.Eye.Position[0] = 0;
            //rt.Eye.Position[1] = 50;
            // big sphere in the back
            Sphere sphBig = rt.AddSphere(0, 100, -600, 100);
            sphBig.AmbientColor = new Vector(0.2, 0.2, 0.2);

            Sphere sphLeftMedium = rt.AddSphere(-60, 75, -200, 50);
            sphLeftMedium.AmbientColor = new Vector(0.2, 0, 0);
            sphLeftMedium.DiffuseColor = new Vector(1.0, 0, 0);

            Sphere sphRightMedium = rt.AddSphere(60, 75, -200, 50);
            sphRightMedium.AmbientColor = new Vector(0, 0.2, 0);

            Sphere sphLeftSmall = rt.AddSphere(-85, 0, -100, 25);
            sphLeftSmall.AmbientColor = new Vector(0, 0, 0.2);

            Sphere sphRightSmall = rt.AddSphere(85, 0, -100, 25);
            sphRightSmall.AmbientColor = new Vector(0, 0.2, 0.2);

            Sphere sphBottom = rt.AddSphere(0, -100, -200, 75);
            sphRightSmall.AmbientColor = new Vector(0.2, 0, 0.2);


            Sphere sphBig2 = rt.AddSphere(0, 100, 800, 100);
            sphBig2.AmbientColor = new Vector(0.2, 0.2, 0.2);

            Sphere sphLeftMedium2 = rt.AddSphere(-60, 75, 400, 50);
            sphLeftMedium2.AmbientColor = new Vector(0.2, 0, 0);
            sphLeftMedium2.DiffuseColor = new Vector(1.0, 0, 0);

            Sphere sphRightMedium2 = rt.AddSphere(60, 75, 400, 50);
            sphRightMedium2.AmbientColor = new Vector(0, 0.2, 0);

            Sphere sphLeftSmall2 = rt.AddSphere(-85, 0, 300, 25);
            sphLeftSmall2.AmbientColor = new Vector(0, 0, 0.2);

            Sphere sphRightSmall2 = rt.AddSphere(85, 0, 300, 25);
            sphRightSmall2.AmbientColor = new Vector(0, 0.2, 0.2);

            //Sphere sphBottom2 = rt.AddSphere(0, -100, 200, 75);
            //sphRightSmall2.AmbientColor = new Vector(0.2, 0, 0.2);


            rt.AddPointLightSource(0, 100, -100);
            rt.AddPointLightSource(200, 200, 0);
            rt.AddPointLightSource(0, 0, -1000);

            pbOutput.Image = rt.Render(pbOutput.Width, pbOutput.Height);
        }


    }
}
