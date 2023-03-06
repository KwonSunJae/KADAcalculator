const express = require("express");
const server = express();
const bodyParser = require('body-parser');
const multer = require('multer');
const form_data = multer();
const fs = require('fs');
const path = require('path');
const spawn = require('child_process').spawn;

server.use(bodyParser.json());
server.use(bodyParser.urlencoded({limit: '5000mb', extended: true, parameterLimit: 100000000000}));

function createoutputhtml (filename,json){
  const ouputfile = path.join(__dirname,path.join('public',filename.toString()+"out.html"));
  html = `
  <!DOCTYPE html>
<html lang="ko">

<head>
    <title>KADA</title>
    <style>
        img{
		float:right;
	}
        div.element {
            width: 410px;
            height: 30px;
        }

        p {
            width: 300px;
            float: left;
            margin: 1px;
        }

        p#fix {
            
            width: 100px;
            float: left;
            
        }
    </style>
</head>

<body>
    <div class="wrapper">
        <form action="./calculate" method="post">

            <div class="title">
                <h1 style="font-size: 21px;">Configuration</h1>
            </div>
            <div class="wrap">
		<img src="${"http://203.252.161.197/public/"+filename.toString()+"trim.png"}"
                <div class="element">
                    <p>fus_apex_FS</p>
                    <p id="fix">${json.fus_apex_FS}</p>
                </div>
                <div class="element">
                    <p>fus_apex_WL</p>
                    <p id="fix">${json.fus_apex_WL}</p>
                </div>
                <div class="element">
                    <p>fus_diameter</p>
                    <p id="fix">${json.fus_diameter}</p>
                </div>
                <div class="element">
                    <p>fus_length</p>
                    <p id="fix">${json.fus_length}</p>
                </div>
                <div class="element">
                    <p>fus_nose_to_length_ratio</p>
                    <p id="fix">${json.fus_nose_to_length_ratio}</p>
                </div>
                <div class="element">
                    <p>fus_tail_to_length_ratio</p>
                    <p id="fix">${json.fus_nose_to_length_ratio}</p>
                </div>

                <div class="element">
                    <p>wing_area</p>
                    <p id="fix">${json.wing_area}</p>
                </div>
                <div class="element">
                    <p>wing_aspect_ratio</p>
                    <p id="fix">${json.wing_aspect_ratio}</p>
                </div>
                <div class="element">
                    <p>wing_taper_ratio</p>
                    <p id="fix">${json.wing_taper_ratio}</p>
                </div>
                <div class="element">
                    <p>wing_leading_edge_sweep</p>
                    <p id="fix">${json.wing_leading_edge_sweep}</p>
                </div>
                <div class="element">
                    <p>wing_apex_FS</p>
                    <p id="fix">${json.wing_apex_FS}</p>
                </div>
                <div class="element">
                    <p>wing_apex_WL</p>
                    <p id="fix">${json.wing_apex_WL}</p>
                </div>
                <div class="element">
                    <p>wing_incidence_angle</p>
                    <p id="fix">${json.wing_incidence_angle}</p>
                </div>
                <div class="element">
                    <p>wing_airfoil_filename</p>
                    <p id="fix">${json.wing_airfoil_filename}</p>
                </div>
                <div class="element">
                    <p>hstab_area</p>
                    <p id="fix">${json.hstab_area}</p>
                </div>
                <div class="element">
                    <p>hstab_aspect_ratio</p>
                    <p id="fix">${json.hstab_aspect_ratio}</p>
                </div>
                <div class="element">
                    <p>hstab_taper_ratio</p>
                    <p id="fix">${json.hstab_taper_ratio}</p>
                </div>
                <div class="element">
                    <p>hstab_leading_edge_sweep</p>
                    <p id="fix">${json.hstab_leading_edge_sweep}</p>
                </div>
                <div class="element">
                    <p>hstab_apex_FS</p>
                    <p id="fix">${json.hstab_apex_FS}</p>
                </div>
                <div class="element">
                    <p>hstab_apex_WL</p>
                    <p id="fix">${json.hstab_apex_WL}</p>
                </div>
                <div class="element">
                    <p>hstab_incidence_angle</p>
                    <p id="fix">${json.hstab_incidence_angle}</p>
                </div>
                <div class="element">
                    <p>hstab_airfoil_filename</p>
                    <p id="fix">${json.hstab_airfoil_filename}</p>
                </div>
                <div class="element">
                    <p>elevator_to_chord_ratio</p>
                    <p id="fix">${json.elevator_to_chord_ratio}</p>
                </div>
                <div class="element">
                    <p>vstab_area</p>
                    <p id="fix">${json.vstab_area}</p>
                </div>
                <div class="element">
                    <p>vstab_aspect_ratio</p>
                    <p id="fix">${json.vstab_aspect_ratio}</p>
                </div>
                <div class="element">
                    <p>vstab_taper_ratio</p>
                    <p id="fix">${json.vstab_taper_ratio}</p>
                </div>
                <div class="element">
                    <p>vstab_leading_edge_sweep</p>
                    <p id="fix">${json.vstab_leading_edge_sweep}</p>
                </div>
                <div class="element">
                    <p>vstab_apex_FS</p>
                    <p id="fix">${json.vstab_apex_FS}</p>
                </div>
                <div class="element">
                    <p>vstab_apex_WL</p>
                    <p id="fix">${json.vstab_apex_WL}</p>
                </div>
                <div class="element">
                    <p>vstab_airfoil_filename</p>
                    <p id="fix">${json.vstab_airfoil_filename}</p>
                </div>

            </div>

            <div class="title">
                <h1 style="font-size: 21px;">Weight</h1>
            </div>
            <div class="wrap">
                <div class="element">
                    <p>mass</p>
                    <p id="fix">${json.mass}</p>
                </div>
                <div class="element">
                    <p>cg_percent_of_mac</p>
                    <p id="fix">${json.cg_percent_of_mac}</p>
                </div>

            </div>
            <div class="title">
                <h1 style="font-size: 21px;">Aerodynamics</h1>
            </div>
            <div class="wrap">
                <div class="element">
                    <p>velocity</p>
                    <p id="fix">${json.velocity}</p>
                </div>
                <div class="element">
                    <p>altitude</p>
                    <p id="fix">${json.altitude}</p>
                </div>
                <div class="element">
                    <p>aoa_min</p>
                    <p id="fix">${json.aoa_min}</p>
                </div>
                <div class="element">
                    <p>aoa_max</p>
                    <p id="fix">${json.aoa_max}</p>
                </div>
                <div class="element">
                    <p>aoa_step</p>
                    <p id="fix">${json.aoa_step}</p>
                </div>

            </div>
            <div class="title">
                <h1 style="font-size: 21px;">Trim</h1>
            </div>
            <div class="wrap">
                <div class="element">
                    <p>velocity_min</p>
                    <p id="fix">${json.velocity_min}</p>
                </div>
                <div class="element">
                    <p>velocity_max</p>
                    <p id="fix">${json.velocity_max}</p>
                </div>
                <div class="element">
                    <p>velocity_step</p>
                    <p id="fix">${json.velocity_step}</p>
                </div>

            </div>
            


            <div class="line">
                <hr>
            </div>

            <div class="submit">
		<a href="../../">HOME </a>
                <button name="home" value="home" onClick="location.href=http://203.252.161.197/">
            </div>
        </form>

    </div>

</body>

</html>
  `;
  fs.writeFile(ouputfile,html,error =>console.log(error));


  return filename+'out.html'

};
const delay = ms => new Promise(resolve => setTimeout(resolve, ms));
server.use('/public',express.static(__dirname + "/public"));

server.get("/", (req, res) => {
 
  res.sendFile(__dirname + "/index.html");
});

server.get("/about", (req, res) => {
  
  res.sendFile(__dirname + "/about.html");
});
server.post("/calculate",async (req,res)=> {
    
    const filename = Date.now();
    let content ='';
    const json = req.body;
    Object.keys(json).forEach(key => {
        content+= key +" = "+ json[key]+"\n";
        
        return key;
    });
    const inputfile = path.join(__dirname,path.join('input',filename+'.in'));
    fs.writeFile(inputfile,content,error =>
        {if (error){
            console.error(error);
        }
    });
    const ouputfile = path.join(__dirname,path.join('public',filename.toString()));

    const result = spawn('python', [path.join(__dirname,'integrated_analysis.py'),inputfile,ouputfile]);
    result.stderr.on('data', function(data) {
        console.log(data.toString());
	res.sendStatus(401);
    });
    await delay(7000);
    return res.redirect('http://203.252.161.197/public/'+createoutputhtml(filename,json));
    
    
});
server.use((req, res) => {
  res.sendFile(__dirname + "/404.html");
});

server.listen(80, (err) => {
  if (err) return console.log(err);
  console.log("The server is listening on port 3000");
});