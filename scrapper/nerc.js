
var webPage = require('webpage');
var page = webPage.create();
var pathnerc = 'scrapper/nerc.html';
var webnerc = 'http://www.nerc.ac.uk/funding/application/currentopportunities/';

var fs = require('fs');
page.open(webnerc, function (status) {
  var content = page.content;
  fs.write(pathnerc, content,'w');
  phantom.exit();
});
