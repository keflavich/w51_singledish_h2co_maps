import keyring
import requests
import bs4

payload = {
    "vanillaLoginForm": "vanillaLoginForm",
    "javax.faces.ViewState": "",
    "ice.window": "",
    "ice.view": "",
    "vanillaLoginForm:vdcId":"",
    "vanillaLoginForm:j_idt23":"",
    "vanillaLoginForm:username":"",
    "vanillaLoginForm:password":"",
    "vanillaLoginForm:affiliateName":"",
    'vanillaLoginForm:button1':'Log in',
    "icefacesCssUpdates": "",
    "vanillaLoginForm:j_idcl":"",
    "javax.faces.source": "vanillaLoginForm:button1",
    "javax.faces.partial.event": "click",
    "javax.faces.partial.execute": "@all",
    "javax.faces.partial.render": "@all",
    "ice.window": "",
    "ice.view": "",
    "ice.focus": "vanillaLoginForm:button1",
    "ice.event.target": "vanillaLoginForm:button1",
    "ice.event.captured": "vanillaLoginForm:button1",
    "ice.event.type": "onclick",
    "ice.event.alt": "false",
    "ice.event.ctrl": "false",
    "ice.event.shift": "false",
    "ice.event.meta": "false",
    "ice.event.left": "true",
    "ice.event.right": "false",
    "ice.submit.type": "ice.s",
    "ice.event.x":963,
    "ice.event.y":295,
    "ice.submit.serialization": "form",
    "javax.faces.partial.ajax": "true",
}

studyid = 99701
S = requests.Session()
r1 = S.get('http://thedata.harvard.edu/dvn/dv/W51_H2CO/faces/login/LoginPage.xhtml')
b1 = bs4.BeautifulSoup(r1.content)

pw = keyring.get_password("http://thedata.harvard.edu/dvn/dv/W51_H2CO/faces/login/LoginPage.xhtml", 'adamginsburg')
payload['javax.faces.ViewState'] = b1.find(id="javax.faces.ViewState").attrs['value']
payload['vanillaLoginForm:vdcId'] = b1.find(id='vanillaLoginForm:vdcId').attrs['value']
payload["vanillaLoginForm:j_idt23"] = b1.find(id="vanillaLoginForm:j_idt23").attrs['value']
affilname = b1.find(id="vanillaLoginForm:affiliateName")
payload["vanillaLoginForm:affiliateName"] = affilname.attrs['value'] if 'value' in affilname.attrs else ''
payload['ice.window'] = b1.find('input', {'name':'ice.window'}).attrs['value']
payload['ice.view'] = b1.find('input', {'name':'ice.view'}).attrs['value']
payload['vanillaLoginForm:username'] = 'adamginsburg'
payload['vanillaLoginForm:password'] = pw

r2 = S.post('http://thedata.harvard.edu/dvn/faces/login/LoginPage.xhtml', data=payload)

r3 = S.get('http://thedata.harvard.edu/dvn/faces/study/TermsOfUsePage.xhtml')


payload2 = {
    "form1": "form1",
    "javax.faces.ViewState": "",
    "ice.window": "",
    "ice.view": "",
    "form1vdcId": "",
    "pageName": "TermsOfUsePage",
    "form1:studyId":studyid,
    "form1:redirectPage":"",
    "form1:tou":"download",
    "form1:termsAccepted":"on",
    "icefacesCssUpdates": "",
    "javax.faces.source": "form1:termsButton",
    "javax.faces.partial.event": "click",
    "javax.faces.partial.execute": "@all",
    "javax.faces.partial.render": "@all",
    "ice.window": "",
    "ice.view": "v",
    "ice.focus": "form1:termsButton",
    "form1": "termsButton:Continue",
    "ice.event.target": "form1:termsButton",
    "ice.event.captured": "form1:termsButton",
    "ice.event.type": "onclick",
    "ice.event.alt": "false",
    "ice.event.ctrl": "false",
    "ice.event.shift": "false",
    "ice.event.meta": "false",
    "ice.event.x": "275",
    "ice.event.y": "688",
    "ice.event.left": "true",
    "ice.event.right": "false",
    "ice.submit.type": "ice.s",
    "ice.submit.serialization": "form",
    "javax.faces.partial.ajax": "true",
}
payload2['javax.faces.ViewState'] = b2.find(id="javax.faces.ViewState").attrs['value']
payload2['form1:vdcId'] = b2.find(id='form1:vdcId').attrs['value']
payload2['ice.window'] = b2.find('input', {'name':'ice.window'}).attrs['value']
payload2['ice.view'] = b2.find('input', {'name':'ice.view'}).attrs['value']

r3 = S.post('http://thedata.harvard.edu/dvn/faces/study/TermsOfUsePage.xhtml', data=payload2)

data = S.get('http://thedata.harvard.edu/dvn/dv/W51_H2CO/faces/study/StudyPage.xhtml?studyId={studyid}&versionNumber=1&checkTermsOfUse=false&cid=6'.format(studyid=studyid))


# VGPS data
"""
http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/data/pub/VGPS/MOS_049_contincluded.Tb?runid=xqby5bqj55qx1uih
http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/data/pub/VGPS/MOS_049.Tb?runid=xqby5bqj55qx1uih
http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/data/pub/VGPS/MOS_049_cont.Tb?runid=xqby5bqj55qx1uih
"""
