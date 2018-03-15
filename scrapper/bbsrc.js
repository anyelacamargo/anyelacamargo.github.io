console.log('Scrapping!');

// scrape_techstars.js

var webPage = require('webpage');
var page = webPage.create();
var fs = require('fs');
var path = 'bbsrc.html';
var webbbsrc = 'https://www.bbsrc.ac.uk/funding/filter/';


page.open(webbbsrc, function (status) {
  var content = page.content;
  fs.write(path, content,'w');
  phantom.exit();
});






