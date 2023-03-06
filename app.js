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
    return res.redirect('http://203.252.161.197/public/'+filename.toString()+'trim.png');
    
    
});
server.use((req, res) => {
  res.sendFile(__dirname + "/404.html");
});

server.listen(80, (err) => {
  if (err) return console.log(err);
  console.log("The server is listening on port 3000");
});