// scrape_techstars.js

var webPage = require('webpage');
var page = webPage.create();

var fs = require('fs');
var path = 'techstars.html';

page.open('https://www.bbsrc.ac.uk/funding/filter/', function (status) {
  var content = page.content;
  fs.write(path,content,'w');
  phantom.exit();
});